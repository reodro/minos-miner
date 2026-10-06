package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadCachingIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * Given a {@link MultiIntervalShard} of {@link GATKRead}, iterates over each {@link AssemblyRegion} within
 * that shard, using the provided {@link AssemblyRegionEvaluator} to determine the boundaries between assembly
 * regions.
 *
 * Loads the reads from the shard as lazily as possible to minimize memory usage.
 *
 * This iterator represents the core of the {@link AssemblyRegionWalker} traversal.
 *
 * NOTE: the provided shard must have appropriate read filters set on it for this traversal type (ie., unmapped
 * and malformed reads must be filtered out).
 */
public class AssemblyRegionIterator implements Iterator<AssemblyRegion> {
    private static final Logger logger = LogManager.getLogger(AssemblyRegionIterator.class);

    private final MultiIntervalShard<GATKRead> readShard;
    private final SAMFileHeader readHeader;
    private final ReferenceDataSource reference;
    private final FeatureManager features;
    private final AssemblyRegionEvaluator evaluator;
    final AssemblyRegionArgumentCollection assemblyRegionArgs;
    
    private AssemblyRegion readyRegion;
    private Queue<AssemblyRegion> pendingRegions;
    private List<GATKRead> previousRegionReads;
    private final ReadCachingIterator readCachingIterator;
    private Queue<GATKRead> readCache;
    private final Iterator<AlignmentContext> locusIterator;     // minos: null when the activity states come from the cache
    private final LocusIteratorByState libs;                   // (ditto)
    private final ActivityProfile activityProfile;
    private Queue<AlignmentAndReferenceContext> pendingAlignmentData;

    /**
     * Constructs an AssemblyRegionIterator over a provided read shard
     *  @param readShard MultiIntervalShard containing the reads that will go into the assembly regions.
     *                  Must have a MAPPED filter set on it.
     * @param readHeader header for the reads
     * @param reference source of reference bases (may be null)
     * @param features source of arbitrary features (may be null)
     * @param evaluator evaluator used to determine whether a locus is active
     */
    public AssemblyRegionIterator(final MultiIntervalShard<GATKRead> readShard,
                                  final SAMFileHeader readHeader,
                                  final ReferenceDataSource reference,
                                  final FeatureManager features,
                                  final AssemblyRegionEvaluator evaluator,
                                  final AssemblyRegionArgumentCollection assemblyRegionArgs,
                                  final boolean trackPileups ) {

        Utils.nonNull(readShard);
        Utils.nonNull(readHeader);
        Utils.nonNull(evaluator);
        assemblyRegionArgs.validate();


        this.readShard = readShard;
        this.readHeader = readHeader;
        this.reference = reference;
        this.features = features;
        this.evaluator = evaluator;
        this.assemblyRegionArgs = assemblyRegionArgs;

        this.readyRegion = null;
        this.previousRegionReads = null;
        this.pendingRegions = new ArrayDeque<>();
        this.readCachingIterator = new ReadCachingIterator(readShard.iterator());
        this.readCache = new ArrayDeque<>();
        this.activityProfile = new BandPassActivityProfile(assemblyRegionArgs.maxProbPropagationDistance, assemblyRegionArgs.activeProbThreshold, BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, readHeader);
        this.pendingAlignmentData = trackPileups ? new ArrayDeque<>() : null;

        final long draws0 = minosDraws();
        final java.nio.file.Path cached = minosActivityFile(trackPileups);
        if ( cached != null && java.nio.file.Files.exists(cached) ) {
            // minos activity cache: this shard's per-locus states were saved by an earlier run of the same skeleton;
            // no pileups, no isActive - the reads are only collected for the regions
            this.libs = null;
            this.locusIterator = null;
            minosReplay = minosOpenReplay(cached);
            minosReads = new htsjdk.samtools.util.PeekableIterator<>(readCachingIterator);
            minosNoDraws(draws0);
        } else {
            // We wrap our LocusIteratorByState inside an IntervalAlignmentContextIterator so that we get empty loci
            // for uncovered locations. This is critical for reproducing GATK 3.x behavior!
            this.libs = new LocusIteratorByState(readCachingIterator, DownsamplingMethod.NONE, ReadUtils.getSamplesFromHeader(readHeader), readHeader, true);
            final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(readShard.getIntervals().iterator());
            this.locusIterator = new IntervalAlignmentContextIterator(libs, intervalLocusIterator, readHeader.getSequenceDictionary());
            if ( cached != null ) {
                minosOpenRecord(cached);
            }
        }
        minosReadDraws += minosDraws() - draws0;

        minosInitStreams(trackPileups);
        readyRegion = loadNextAssemblyRegion();
    }

    @Override
    public boolean hasNext() {
        return readyRegion != null;
    }

    @Override
    public AssemblyRegion next() {
        if (!hasNext()) {
            throw new NoSuchElementException("next() called when there were no more elements");
        }

        final AssemblyRegion toReturn = readyRegion;
        previousRegionReads = toReturn.getReads();
        readyRegion = loadNextAssemblyRegion();
        return toReturn;
    }

    private AssemblyRegion loadNextAssemblyRegion() {
        AssemblyRegion nextRegion = null;
        final long draws0 = minosDraws();
        final long excluded0 = minosExcludedDraws;

        // minos: the loop condition is evaluated exactly as before (hasNext first), and variant streams complete a
        // region found at the previous locus right after it - the point where their own run would return it.
        while ( true ) {
            final boolean moreLoci = minosMoreLoci();
            if ( ! minosStreams.isEmpty() ) {
                final long d = minosDraws();             // variant streams genotype on their own generators
                minosAdvanceStreams(moreLoci);
                minosExcludedDraws += minosDraws() - d;
            }
            if ( ! (moreLoci && nextRegion == null) ) {
                break;
            }
            // minos activity cache: the locus and its state come from the cache, or from the pileup and isActive
            final AlignmentContext pileup;
            final SimpleInterval pileupInterval;
            ActivityProfileState state = null;
            if ( minosReplay != null ) {
                state = minosReplayState();
                pileup = null;
                pileupInterval = state.getLoc();
            } else {
                pileup = locusIterator.next();
                logger.debug(()->"AssemblyRegionIterator pileup: " + pileup);
                pileupInterval = new SimpleInterval(pileup);
            }

            // Pop any new pending regions off of the activity profile. These pending regions will not become ready
            // until we've traversed all the reads that belong in them.
            //
            // Ordering matters here: need to check for forceConversion before adding current pileup to the activity profile
            if ( ! activityProfile.isEmpty() ) {
                final boolean forceConversion = pileupInterval.getStart() != activityProfile.getEnd() + 1;
                pendingRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionArgs.assemblyRegionPadding, assemblyRegionArgs.minAssemblyRegionSize, assemblyRegionArgs.maxAssemblyRegionSize, forceConversion));
            }

            // Add the current pileup to the activity profile
            ReferenceContext pileupRefContext = null;
            if ( state == null ) {
                pileupRefContext = new ReferenceContext(reference, pileupInterval);
                final FeatureContext pileupFeatureContext = new FeatureContext(features, pileupInterval);
                state = evaluator.isActive(pileup, pileupRefContext, pileupFeatureContext);
                minosRecordState(state);                 // before the profile stores (and later modifies) it
            }
            final ActivityProfileState profile = state;
            logger.debug(() -> profile.toString());
            // minos: profiles add neighbours' smoothed probability onto the state objects they store, so every
            // stream gets its own copy, taken before the main profile stores (and later modifies) the original
            final List<ActivityProfileState> minosStates = new ArrayList<>(minosStreams.size());
            for ( int i = 0; i < minosStreams.size(); i++ ) {
                minosStates.add(minosCopyState(profile));
            }
            activityProfile.add(profile);

            if (pendingAlignmentData!=null) {
                AlignmentAndReferenceContext alignmentAndRefContext = new AlignmentAndReferenceContext(pileup, pileupRefContext);
                alignmentAndRefContext.setActivityScore(profile.getOriginalActiveProb());
                pendingAlignmentData.add(alignmentAndRefContext);
            }

            // A pending region only becomes ready once our locus iterator has advanced beyond the end of its extended span
            // (this ensures that we've loaded all reads that belong in the new region)
            if ( ! pendingRegions.isEmpty() && IntervalUtils.isAfter(pileupInterval, pendingRegions.peek().getPaddedSpan(), readHeader.getSequenceDictionary()) ) {
                nextRegion = pendingRegions.poll();
            }
            for ( int i = 0; i < minosStreams.size(); i++ ) {
                minosStreams.get(i).onLocus(pileupInterval, minosStates.get(i));
            }
        }

        // When we run out of loci, close out the activity profile, and close out any remaining pending regions one at a time
        // It may require multiple invocations before the pendingRegions queue is cleared out.
        if ( ! minosMoreLoci() ) {

            // Pull on the encapsulated LocusIteratorByState until it's exhausted, as the enclosing
            // IntervalLocusIterator is not guaranteed to exhaust it. This guarantees that the reads
            // in the final padded region end up in our read cache.
            minosDrainReads();
            minosFinishShard((minosDraws() - draws0) - (minosExcludedDraws - excluded0));

            if ( ! activityProfile.isEmpty() ) {
                // Pop the activity profile a final time with forceConversion == true
                pendingRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionArgs.assemblyRegionPadding, assemblyRegionArgs.minAssemblyRegionSize, assemblyRegionArgs.maxAssemblyRegionSize, true));
            }

            // Grab the next pending region if there is one, unless we already have a region ready to go
            // (could happen if we finalize a region on the same iteration that we run out of loci in the locus iterator)
            if ( ! pendingRegions.isEmpty() && nextRegion == null ) {
                nextRegion = pendingRegions.poll();
            }
        }

        // If there's a region ready, fill it with reads before returning
        if ( nextRegion != null ) {
            fillNextAssemblyRegionWithReads(nextRegion);
            // fillnextessemblyregion; check you are on correct chr; if alignment data is not in the assembly region then pop it
            fillNextAssemblyRegionWithPileupData(nextRegion);
        }

        minosReadDraws += (minosDraws() - draws0) - (minosExcludedDraws - excluded0);
        return nextRegion;
    }

    private void fillNextAssemblyRegionWithReads( final AssemblyRegion region ) {
        // First we need to check the previous region for reads that also belong in this region
        if ( previousRegionReads != null ) {
            for ( final GATKRead previousRegionRead : previousRegionReads ) {
                if ( region.getPaddedSpan().overlaps(previousRegionRead) ) {
                    region.add(previousRegionRead);
                }
            }
        }

        // Update our read cache with the reads newly-read by our locus iterator.
        // Note that the cache is (crucially) maintained in coordinate order.
        minosConsumeReads();

        // Add all reads from the cache that belong in this region
        while ( ! readCache.isEmpty() ) {
            final GATKRead nextRead = readCache.peek();

            // Stop once we encounter a read that starts after the end of the last region's span
            // (and leave it in the readCache)
            if ( IntervalUtils.isAfter(nextRead, region.getPaddedSpan(), readHeader.getSequenceDictionary()) ) {
                break;
            }

            // Ok, we're going to consume this read
            readCache.poll();

            // Add the read if it overlaps the region's extended span. If it doesn't, it must end before the
            // start of the region's extended span, so we discard it.
            if ( region.getPaddedSpan().overlaps(nextRead) ) {
                region.add(nextRead);
            }
        }
    }

    private void fillNextAssemblyRegionWithPileupData(final AssemblyRegion region){
        // Save ourselves the memory footprint and work of saving the pileups in the event they aren't needed for processing.
        if (pendingAlignmentData == null){
            return;
        }
        final List<AlignmentAndReferenceContext> overlappingAlignmentData = new ArrayList<>();
        final Queue<AlignmentAndReferenceContext> previousAlignmentData = new ArrayDeque<>();

        while (!pendingAlignmentData.isEmpty()) {
            final AlignmentContext pendingAlignmentContext = pendingAlignmentData.peek().getAlignmentContext();
            if (!pendingAlignmentContext.contigsMatch(region) ||
                pendingAlignmentContext.getStart() < region.getStart()) {
                pendingAlignmentData.poll(); // pop this
            } else {
                break;
            }
        }
        while (!pendingAlignmentData.isEmpty()) {
            final AlignmentContext pendingAlignmentContext = pendingAlignmentData.peek().getAlignmentContext();

            if (!pendingAlignmentContext.contigsMatch(region) ||
                    pendingAlignmentContext.getStart() <= region.getEnd()) {
                overlappingAlignmentData.add(pendingAlignmentData.poll()); // pop into overlappingAlignmentData
            } else {
                break;
            }
        }

        // reconstructing queue to contain items that may be in the next assembly region
        previousAlignmentData.addAll(overlappingAlignmentData);
        previousAlignmentData.addAll(pendingAlignmentData);
        pendingAlignmentData = previousAlignmentData;

        region.addAllAlignmentData(overlappingAlignmentData);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported by AssemblyRegionIterator");
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos: activity-profile branching. Config variants that differ from the main config only in
    // active_probability_threshold / min_assembly_region_size / max_assembly_region_size see the same loci, the same
    // per-locus activity states (isActive does not read these settings) and the same reads. Each gets its own
    // activity profile, pending regions and read cache fed from the shared traversal, and its regions go to the
    // engine at the moment its own run would process them: right after the stream's next region is loaded.
    // ---------------------------------------------------------------------------------------------------------------
    public interface MinosActivityHook {
        /** Activity settings {threshold, min size, max size} of each variant stream that differs from the main one. */
        List<double[]> minosActivityStreams(double activeProbThreshold, int minRegionSize, int maxRegionSize);

        /** Call variants for one region of a stream, as HaplotypeCaller.apply would in that variant's own run. */
        void minosProcessRegion(int stream, AssemblyRegion region, ReferenceContext referenceContext, FeatureContext featureContext);
    }

    public static volatile MinosActivityHook minosHook = null;

    private final List<MinosStream> minosStreams = new ArrayList<>();

    private void minosInitStreams(final boolean trackPileups) {
        final MinosActivityHook hook = minosHook;
        if ( hook == null ) {
            return;
        }
        final List<double[]> specs = hook.minosActivityStreams(assemblyRegionArgs.activeProbThreshold,
                assemblyRegionArgs.minAssemblyRegionSize, assemblyRegionArgs.maxAssemblyRegionSize);
        if ( ! specs.isEmpty() && (trackPileups || assemblyRegionArgs.forceActive) ) {
            throw new IllegalStateException("minos: activity variants do not support pileup tracking or forced activity");
        }
        for ( int i = 0; i < specs.size(); i++ ) {
            minosStreams.add(new MinosStream(i, hook, specs.get(i)[0], (int) specs.get(i)[1], (int) specs.get(i)[2]));
        }
    }

    private static ActivityProfileState minosCopyState(final ActivityProfileState state) {
        final ActivityProfileState copy = new ActivityProfileState(state.getLoc(), state.isActiveProb(), state.getResultState(), state.getResultValue());
        copy.setOriginalActiveProb(state.getOriginalActiveProb());
        return copy;
    }

    /** Reads newly loaded by the locus iterator go, in the same coordinate order, to every read cache. */
    private void minosConsumeReads() {
        final List<GATKRead> newReads = readCachingIterator.consumeCachedReads();
        readCache.addAll(newReads);
        for ( final MinosStream s : minosStreams ) {
            s.readCache.addAll(newReads);
        }
    }

    private void minosAdvanceStreams(final boolean moreLoci) {
        if ( ! moreLoci ) {
            minosDrainReads();               // the same drain the main path does when loci run out; idempotent
        }
        for ( final MinosStream s : minosStreams ) {
            s.advance(moreLoci);
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos: activity cache (-Dminos.activity.dir=<one directory per skeleton>). The per-locus activity states depend
    // only on what a batch shares (the skeleton: everything but the branchable parameters), so a run saves each shard's
    // states and a later run of the same skeleton reads them instead of building pileups and calling isActive (about
    // half of a run). The reads are still collected, now ahead of the loci instead of with the pileups: every read
    // starting at or before a locus is in the cache when that locus is processed, so each region gets the same reads.
    // Collecting them at other moments is exact only if the read pass draws no random numbers (positional downsampling
    // at deep spots shares java's Random with the QD jitter of genotyping), so a shard is saved only when its read pass
    // drew nothing - and a read from the cache that would draw fails the run instead of returning different calls.
    // ---------------------------------------------------------------------------------------------------------------
    private static final int MINOS_ACT_MAGIC = 0x4d414354;        // "MACT"
    private static final int MINOS_ACT_VERSION = 1;
    private static final ActivityProfileState.Type[] MINOS_ACT_TYPES = ActivityProfileState.Type.values();

    private DataInputStream minosReplay = null;                   // reading: the shard's saved states
    private ActivityProfileState minosNext = null;                // ... the next one (null: no more loci)
    private htsjdk.samtools.util.PeekableIterator<GATKRead> minosReads = null;
    private String minosContig = null;
    private DataOutputStream minosOut = null;                     // saving: the states of this run
    private java.nio.file.Path minosOutTmp = null;
    private java.nio.file.Path minosOutPath = null;
    private boolean minosFinished = false;
    private long minosReadDraws = 0;                              // raw draws of java's Random during the read pass
    private long minosExcludedDraws = 0;                          // ... minus the variant streams' own genotyping

    private static long minosDraws() {
        return Utils.minosGlobalRandom().minosRawDraws();
    }

    private java.nio.file.Path minosActivityFile(final boolean trackPileups) {
        final String dir = System.getProperty("minos.activity.dir");
        if ( dir == null || dir.isEmpty() || trackPileups || assemblyRegionArgs.forceActive || readShard.getIntervals().size() != 1 ) {
            return null;
        }
        final SimpleInterval iv = readShard.getIntervals().get(0);
        return java.nio.file.Paths.get(dir, iv.getContig() + "_" + iv.getStart() + "_" + iv.getEnd() + ".act.gz");
    }

    private boolean minosMoreLoci() {
        return minosReplay != null ? minosNext != null : locusIterator.hasNext();
    }

    private void minosDrainReads() {
        if ( minosReplay == null ) {
            while ( libs.hasNext() ) {
                libs.next();
            }
            return;
        }
        final long d = minosDraws();
        while ( minosReads.hasNext() ) {
            minosReads.next();
        }
        minosNoDraws(d);
    }

    private void minosNoDraws(final long before) {
        if ( minosDraws() != before ) {
            throw new IllegalStateException("minos activity cache: the read pass drew random numbers on a cached shard");
        }
    }

    private void minosOpenRecord(final java.nio.file.Path path) {
        try {
            java.nio.file.Files.createDirectories(path.getParent());
            minosOutPath = path;
            minosOutTmp = path.resolveSibling(path.getFileName() + "." + java.util.UUID.randomUUID() + ".tmp");
            minosOut = new DataOutputStream(new java.io.BufferedOutputStream(
                    new java.util.zip.GZIPOutputStream(java.nio.file.Files.newOutputStream(minosOutTmp), 1 << 16) {
                        { def.setLevel(java.util.zip.Deflater.BEST_SPEED); }
                    }, 1 << 16));
            final SimpleInterval iv = readShard.getIntervals().get(0);
            minosOut.writeInt(MINOS_ACT_MAGIC);
            minosOut.writeInt(MINOS_ACT_VERSION);
            minosOut.writeUTF(iv.getContig());
            minosOut.writeInt(iv.getStart());
            minosOut.writeInt(iv.getEnd());
            minosContig = iv.getContig();
        } catch ( final IOException e ) {
            logger.warn("minos activity cache: cannot save " + path + ": " + e);
            minosAbortRecord();
        }
    }

    private void minosRecordState(final ActivityProfileState s) {
        if ( minosOut == null ) {
            return;
        }
        try {
            if ( ! s.getLoc().getContig().equals(minosContig) || s.getLoc().getStart() != s.getLoc().getEnd() ) {
                minosAbortRecord();
                return;
            }
            final Number value = s.getResultValue();
            minosOut.writeInt(s.getLoc().getStart());
            minosOut.writeDouble(s.isActiveProb());
            minosOut.writeDouble(s.getOriginalActiveProb());
            minosOut.writeByte(s.getResultState().ordinal() | (value != null ? 0x10 : 0));
            if ( value != null ) {
                minosOut.writeDouble(value.doubleValue());
            }
        } catch ( final IOException e ) {
            logger.warn("minos activity cache: saving failed: " + e);
            minosAbortRecord();
        }
    }

    private void minosAbortRecord() {
        try {
            if ( minosOut != null ) {
                minosOut.close();
            }
        } catch ( final IOException ignored ) {
            // the file is deleted below
        }
        minosOut = null;
        if ( minosOutTmp != null ) {
            try {
                java.nio.file.Files.deleteIfExists(minosOutTmp);
            } catch ( final IOException ignored ) {
                // a leftover .tmp is never read
            }
        }
    }

    /** Out of loci (called again while pending regions drain: only the first call acts). */
    private void minosFinishShard(final long drawsInProgress) {
        if ( minosFinished ) {
            return;
        }
        minosFinished = true;
        if ( minosReplay != null ) {
            try {
                minosReplay.close();
            } catch ( final IOException ignored ) {
                // read completely
            }
            return;
        }
        if ( minosOut == null ) {
            return;
        }
        if ( minosReadDraws + drawsInProgress != 0 ) {           // deep spots: this shard's states are not reusable
            minosAbortRecord();
            return;
        }
        try {
            minosOut.writeInt(-1);
            minosOut.close();
            minosOut = null;
            java.nio.file.Files.move(minosOutTmp, minosOutPath, java.nio.file.StandardCopyOption.ATOMIC_MOVE,
                    java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        } catch ( final IOException e ) {
            logger.warn("minos activity cache: cannot save " + minosOutPath + ": " + e);
            minosAbortRecord();
        }
    }

    private DataInputStream minosOpenReplay(final java.nio.file.Path path) {
        try {
            final DataInputStream in = new DataInputStream(new java.io.BufferedInputStream(
                    new java.util.zip.GZIPInputStream(java.nio.file.Files.newInputStream(path), 1 << 16), 1 << 16));
            final SimpleInterval iv = readShard.getIntervals().get(0);
            if ( in.readInt() != MINOS_ACT_MAGIC || in.readInt() != MINOS_ACT_VERSION || ! in.readUTF().equals(iv.getContig())
                    || in.readInt() != iv.getStart() || in.readInt() != iv.getEnd() ) {
                throw new IllegalStateException("minos activity cache: " + path + " is not this shard's");
            }
            minosContig = iv.getContig();
            minosNext = minosReadState(in);
            return in;
        } catch ( final IOException e ) {
            throw new IllegalStateException("minos activity cache: cannot read " + path, e);
        }
    }

    private ActivityProfileState minosReadState(final DataInputStream in) throws IOException {
        final int pos = in.readInt();
        if ( pos < 0 ) {
            return null;
        }
        final double prob = in.readDouble();
        final double original = in.readDouble();
        final int flags = in.readByte();
        final SimpleInterval loc = new SimpleInterval(minosContig, pos, pos);
        final ActivityProfileState.Type type = MINOS_ACT_TYPES[flags & 0x0f];
        final ActivityProfileState s;
        if ( (flags & 0x10) != 0 ) {
            s = new ActivityProfileState(loc, prob, type, in.readDouble());
        } else if ( type == ActivityProfileState.Type.NONE ) {
            s = new ActivityProfileState(loc, prob);             // the two-argument form isActive uses with no data
        } else {
            s = new ActivityProfileState(loc, prob, type, null);
        }
        s.setOriginalActiveProb(original);
        return s;
    }

    /** The next locus's saved state, with every read starting at or before it pulled into the read cache. */
    private ActivityProfileState minosReplayState() {
        final ActivityProfileState s = minosNext;
        final long d = minosDraws();
        while ( minosReads.hasNext() && ! IntervalUtils.isAfter(minosReads.peek(), s.getLoc(), readHeader.getSequenceDictionary()) ) {
            minosReads.next();
        }
        minosNoDraws(d);
        try {
            minosNext = minosReadState(minosReplay);
        } catch ( final IOException e ) {
            throw new IllegalStateException("minos activity cache: cannot read the saved states", e);
        }
        return s;
    }

    private final class MinosStream {
        final int id;
        final MinosActivityHook hook;
        final int minRegionSize;
        final int maxRegionSize;
        final ActivityProfile profile;
        final Queue<AssemblyRegion> pending = new ArrayDeque<>();
        final Queue<GATKRead> readCache = new ArrayDeque<>();
        List<GATKRead> previousRegionReads = null;
        AssemblyRegion ready = null;        // returned by its own next(); processed once the following one is loaded
        AssemblyRegion loaded = null;       // became ready at the last locus, not yet filled

        MinosStream(final int id, final MinosActivityHook hook, final double threshold, final int minRegionSize, final int maxRegionSize) {
            this.id = id;
            this.hook = hook;
            this.minRegionSize = minRegionSize;
            this.maxRegionSize = maxRegionSize;
            this.profile = new BandPassActivityProfile(assemblyRegionArgs.maxProbPropagationDistance, threshold,
                    BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, readHeader);
        }

        /** The main loop body for this stream, on the shared locus and activity state. */
        void onLocus(final SimpleInterval locus, final ActivityProfileState state) {
            if ( ! profile.isEmpty() ) {
                final boolean forceConversion = locus.getStart() != profile.getEnd() + 1;
                pending.addAll(profile.popReadyAssemblyRegions(assemblyRegionArgs.assemblyRegionPadding, minRegionSize, maxRegionSize, forceConversion));
            }
            profile.add(state);
            if ( ! pending.isEmpty() && IntervalUtils.isAfter(locus, pending.peek().getPaddedSpan(), readHeader.getSequenceDictionary()) ) {
                loaded = pending.poll();
            }
        }

        void advance(final boolean moreLoci) {
            if ( moreLoci ) {
                if ( loaded != null ) {
                    complete();
                }
                return;
            }
            // out of loci: final forced pop, then the remaining regions one by one, then the last ready one
            if ( ! profile.isEmpty() ) {
                pending.addAll(profile.popReadyAssemblyRegions(assemblyRegionArgs.assemblyRegionPadding, minRegionSize, maxRegionSize, true));
            }
            if ( loaded == null ) {
                loaded = pending.poll();
            }
            while ( loaded != null ) {
                complete();
                loaded = pending.poll();
            }
            if ( ready != null ) {
                final AssemblyRegion last = ready;
                ready = null;
                process(last);
            }
        }

        /** Own run: next() returned `ready`, then loaded `loaded` (with ready's reads as the previous ones), then apply(ready). */
        private void complete() {
            if ( ready != null ) {
                previousRegionReads = ready.getReads();
            }
            fill(loaded);
            final AssemblyRegion toProcess = ready;
            ready = loaded;
            loaded = null;
            if ( toProcess != null ) {
                process(toProcess);
            }
        }

        private void process(final AssemblyRegion region) {
            hook.minosProcessRegion(id, region, new ReferenceContext(reference, region.getPaddedSpan()),
                    new FeatureContext(features, region.getPaddedSpan()));
        }

        private void fill(final AssemblyRegion region) {
            if ( previousRegionReads != null ) {
                for ( final GATKRead previousRegionRead : previousRegionReads ) {
                    if ( region.getPaddedSpan().overlaps(previousRegionRead) ) {
                        region.add(previousRegionRead);
                    }
                }
            }
            minosConsumeReads();
            while ( ! readCache.isEmpty() ) {
                final GATKRead nextRead = readCache.peek();
                if ( IntervalUtils.isAfter(nextRead, region.getPaddedSpan(), readHeader.getSequenceDictionary()) ) {
                    break;
                }
                readCache.poll();
                if ( region.getPaddedSpan().overlaps(nextRead) ) {
                    region.add(nextRead);
                }
            }
        }
    }
}
