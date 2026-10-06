package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.MinosKBest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.List;
import java.util.UUID;
import java.util.function.Supplier;

/**
 * minos: assembly cache (-Dminos.assembly.dir=<directory of the round>, -Dminos.assembly.skeleton=<the config's
 * skeleton id: every setting outside the branchable stages, e.g. min_base_quality_score, which the assembler reads>).
 *
 * A region's local assembly (AssemblyBasedCallerUtils.assembleReads) depends on the region (span, padding), its
 * finalized reads, the reference and the assembly settings - never on the PairHMM or genotyping settings. Configs of
 * a round that differ only in those later stages assemble the same regions again; here the first assembly of a
 * region saves the haplotypes it found (in order, every field the later steps read), and a later run with the same
 * key rebuilds the result on GATK's own empty assembly result (what runLocalAssembly starts from) with its own
 * reference haplotype and the saved haplotypes in their saved order (replaceAllHaplotypes, as the assembler itself
 * finishes), instead of building and searching the graphs.
 *
 * Key: the skeleton id plus the branchable assembly settings (min pruning, dangling branches, haplotype cap, pruning
 * error rate and LOD), the region's span and padded span, and every finalized read as its SAM text (all
 * fields and tags) in order - so a hit means the assembler would have seen exactly the same input. The k-best searches
 * an assembly reports to MinosKBest (the learned haplotype cap) are saved with it and reported again on a hit.
 * Nothing is saved for flow-based collapsing or partially determined haplotypes (never on this pipeline).
 */
public final class MinosAssemblyCache {
    private static final Logger logger = LogManager.getLogger(MinosAssemblyCache.class);
    private static final int MAGIC = 0x4d415353;        // "MASS"
    private static final int VERSION = 1;

    private MinosAssemblyCache() { }

    /**
     * The assembly of `region` under `assemblyKey` (the settings `full` assembles with): read from the cache, or
     * assembled by `full` and saved. `bypass` is assembleReads with bypassAssembly=true: it finalizes the region
     * (assembleReads then skips that) and returns the reference-only result a hit is built on.
     */
    public static AssemblyResultSet assemble(final String assemblyKey, final AssemblyRegion region, final SAMFileHeader header,
                                             final Supplier<AssemblyResultSet> bypass, final Supplier<AssemblyResultSet> full) {
        final String dir = System.getProperty("minos.assembly.dir");
        final String skeleton = System.getProperty("minos.assembly.skeleton");    // every other setting assembly reads
        if (dir == null || dir.isEmpty() || skeleton == null || skeleton.isEmpty()) {
            return full.get();
        }
        final AssemblyResultSet empty = bypass.get();
        final Path path;
        try {
            final String key = key(skeleton + "|" + assemblyKey, region, header);
            path = Paths.get(dir, key.substring(0, 2), key + ".hap.gz");
        } catch (final Exception e) {
            logger.warn("minos assembly cache: no key for " + region.getSpan() + ": " + e);
            return full.get();
        }
        if (Files.exists(path)) {
            final AssemblyResultSet loaded = load(path, empty);
            if (loaded != null) {
                return loaded;
            }
        }
        MinosKBest.startCapture();
        final AssemblyResultSet result;
        final List<int[]> searches;
        try {
            result = full.get();
        } finally {
            searches = MinosKBest.stopCapture();
        }
        save(path, result, searches);
        return result;
    }

    static String key(final String assemblyKey, final AssemblyRegion region, final SAMFileHeader header) throws Exception {
        final MessageDigest md = MessageDigest.getInstance("SHA-256");
        md.update((VERSION + "|" + assemblyKey + "|" + region.getSpan() + "|" + region.getPaddedSpan() + "|"
                + region.isActive() + "|" + region.getReads().size() + "\n").getBytes(StandardCharsets.UTF_8));
        for (final GATKRead read : region.getReads()) {
            md.update(read.convertToSAMRecord(header).getSAMString().getBytes(StandardCharsets.UTF_8));
            md.update((byte) 0);
        }
        final StringBuilder hex = new StringBuilder();
        for (final byte b : md.digest()) {
            hex.append(String.format("%02x", b));
        }
        return hex.substring(0, 40);
    }

    private static void save(final Path path, final AssemblyResultSet result, final List<int[]> searches) {
        if (result.getHaplotypeCollapsingEngine() != null || result.isPartiallyDeterminedList()) {
            return;
        }
        final List<Haplotype> haplotypes = result.getHaplotypeList();
        for (final Haplotype h : haplotypes) {
            if (h.getClass() != Haplotype.class || h.getGenomeLocation() == null) {
                return;
            }
        }
        Path tmp = null;
        try {
            Files.createDirectories(path.getParent());
            tmp = path.resolveSibling(path.getFileName() + "." + UUID.randomUUID() + ".tmp");
            try (DataOutputStream out = new DataOutputStream(new java.io.BufferedOutputStream(
                    new java.util.zip.GZIPOutputStream(Files.newOutputStream(tmp))))) {
                out.writeInt(MAGIC);
                out.writeInt(VERSION);
                out.writeInt(searches.size());
                for (final int[] s : searches) {
                    out.writeInt(s[0]);
                    out.writeBoolean(s[1] != 0);
                }
                out.writeInt(haplotypes.size());
                for (final Haplotype h : haplotypes) {
                    final byte[] bases = h.getBases();
                    out.writeInt(bases.length);
                    out.write(bases);
                    out.writeBoolean(h.isReference());
                    out.writeUTF(h.getGenomeLocation().getContig());
                    out.writeInt(h.getGenomeLocation().getStart());
                    out.writeInt(h.getGenomeLocation().getEnd());
                    final Cigar cigar = h.getCigar();
                    out.writeBoolean(cigar != null);
                    if (cigar != null) {
                        out.writeUTF(cigar.toString());
                    }
                    out.writeInt(h.getAlignmentStartHapwrtRef());
                    out.writeLong(Double.doubleToRawLongBits(h.getScore()));
                    out.writeInt(h.getKmerSize());
                    out.writeBoolean(h.isCollapsed());
                    out.writeInt(h.getUniquenessValue());
                }
            }
            Files.move(tmp, path, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (final IOException e) {
            logger.warn("minos assembly cache: cannot save " + path + ": " + e);
            if (tmp != null) {
                try {
                    Files.deleteIfExists(tmp);
                } catch (final IOException ignored) {
                    // a leftover .tmp is never read
                }
            }
        }
    }

    /** The saved assembly on top of `empty` (null: not usable, assemble instead). */
    private static AssemblyResultSet load(final Path path, final AssemblyResultSet empty) {
        final List<Haplotype> start = empty.getHaplotypeList();
        if (start.size() != 1 || !start.get(0).isReference()) {
            return null;
        }
        final Haplotype ref = start.get(0);
        try (DataInputStream in = new DataInputStream(new java.io.BufferedInputStream(
                new java.util.zip.GZIPInputStream(Files.newInputStream(path))))) {
            if (in.readInt() != MAGIC || in.readInt() != VERSION) {
                return null;
            }
            final int nSearches = in.readInt();
            final int[][] searches = new int[nSearches][];
            for (int i = 0; i < nSearches; i++) {
                searches[i] = new int[]{in.readInt(), in.readBoolean() ? 1 : 0};
            }
            final int n = in.readInt();
            final Haplotype[] haplotypes = new Haplotype[n];
            boolean refSeen = false;
            for (int i = 0; i < n; i++) {
                final byte[] bases = new byte[in.readInt()];
                in.readFully(bases);
                final boolean isRef = in.readBoolean();
                final SimpleInterval loc = new SimpleInterval(in.readUTF(), in.readInt(), in.readInt());
                final Cigar cigar = in.readBoolean() ? TextCigarCodec.decode(in.readUTF()) : null;
                final Haplotype h;
                if (isRef) {                                  // the reference haplotype (anywhere): GATK's own object
                    if (refSeen || !Arrays.equals(bases, ref.getBases()) || !loc.equals(new SimpleInterval(ref.getGenomeLocation()))) {
                        logger.warn("minos assembly cache: " + path + " has another reference haplotype; assembling");
                        return null;
                    }
                    refSeen = true;
                    h = ref;
                } else {
                    h = new Haplotype(bases, false);
                    h.setGenomeLocation(loc);
                }
                if (cigar != null) {
                    h.setCigar(cigar);
                }
                h.setAlignmentStartHapwrtRef(in.readInt());
                h.setScore(Double.longBitsToDouble(in.readLong()));
                h.setKmerSize(in.readInt());
                h.setCollapsed(in.readBoolean());
                h.setUniquenessValue(in.readInt());
                haplotypes[i] = h;
            }
            if (!refSeen) {
                return null;
            }
            // the saved order exactly (the assembler may end with replaceAllHaplotypes, which puts the reference last and
            // can keep a non-reference haplotype with the reference's bases): rebuilt the way GATK builds such a list
            empty.replaceAllHaplotypes(Arrays.asList(haplotypes));
            if (empty.getHaplotypeList().size() != n) {
                return null;
            }
            for (final int[] s : searches) {                  // the cap report, as if the searches ran here
                MinosKBest.record(s[0], s[1] != 0);
            }
            return empty;
        } catch (final IOException | RuntimeException e) {
            logger.warn("minos assembly cache: cannot read " + path + ": " + e + "; assembling");
            return null;
        }
    }
}
