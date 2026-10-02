package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.NamedFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.OutputMode;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.transformers.IUPACReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * The core engine for the HaplotypeCaller that does all of the actual work of the tool.
 * Usage:
 * -Pass the HaplotypeCaller args into the constructor, which will initialize the HC engine completely.
 * -Get the appropriate VCF or GVCF writer (depending on our arguments) from {@link #makeVCFWriter}
 * -Write the appropriate VCF header via {@link #writeHeader}
 * -Repeatedly call {@link #isActive} to identify active vs. inactive regions
 * -Repeatedly call {@link #callRegion} to call variants in each region, and add them to your writer
 * -When done, call {@link #shutdown}. Close the writer you got from {@link #makeVCFWriter} yourself.
 */
public class HaplotypeCallerEngine implements AssemblyRegionEvaluator, AssemblyRegionIterator.MinosActivityHook {

    private static final Logger logger = LogManager.getLogger(HaplotypeCallerEngine.class);

    protected final HaplotypeCallerArgumentCollection hcArgs;

    protected final SAMFileHeader readsHeader;

    protected ReferenceConfidenceModel referenceConfidenceModel = null;

    protected AssemblyRegionTrimmer trimmer;

    protected final OutputStreamWriter assemblyDebugOutStream;

    /**
     * List of interval file entries including regions and custom ploidy values to apply in that region.
     */
    protected final List<SimpleCount> ploidyRegions = new ArrayList<>();

    /**
     * An OverlapDetector object for checking whether region overlaps given ploidyRegions.
     */
    protected final OverlapDetector<SimpleCount> ploidyRegionsOverlapDetector;

    /**
     * List of all custom ploidies provided by user
     */
    private final LinkedHashSet<Integer> allCustomPloidies;

    /**
     * The default genotyping engine for the isActive() determination
     */
    private MinimalGenotypingEngine defaultActiveRegionEvaluationGenotyperEngine = null;

    /**
     * Map of user-provided ploidy values to corresponding active region genotyper. Values are checked as valid Integers during
     * initialization, but Strings are used as keys to avoid parsing repeatedly during runtime.
     */
    private final HashMap<Integer, MinimalGenotypingEngine> ploidyToActiveEvaluationGenotyper = new HashMap<>();

    protected ReadThreadingAssembler assemblyEngine = null;

    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;
    private ReadLikelihoodCalculationEngine filterStepLikelihoodCalculationEngine = null;
    // If we are in PDHMM mode we need to hold onto two likelihoods engines for the fallback
    private ReadLikelihoodCalculationEngine pdhmmLikelihoodCalculationEngine = null;

    /**
     * The default genotyping engine to use for actual variant calling and genotyping in an active region.
     */
    protected HaplotypeCallerGenotypingEngine defaultGenotypingEngine = null;

    /**
     * Map of user-provided ploidy values to corresponding genotyper. Values are checked as valid Integers during
     * initialization, but Strings are used as keys to avoid parsing repeatedly during runtime.
     */
    protected final HashMap<Integer, HaplotypeCallerGenotypingEngine> ploidyToGenotyperMap = new HashMap<>();

    private VariantAnnotatorEngine annotationEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    protected final ReferenceSequenceFile referenceReader;

    // writes Haplotypes to a bam file when the -bamout option is specified
    protected Optional<HaplotypeBAMWriter> haplotypeBAMWriter;

    // writes Variants from assembly graph
    private Optional<VariantContextWriter> assembledEventMapVcfOutputWriter;
    private Optional<PriorityQueue<VariantContext>> assembledEventMapVariants;

    protected Optional<AlleleLikelihoodWriter> alleleLikelihoodWriter;

    private Set<String> sampleSet;
    protected SampleList samplesList;

    private final boolean forceCallingAllelesPresent;

    private byte minTailQuality;

    protected SmithWatermanAligner aligner;

    private final DragstrParams dragstrParams;

    public static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    private static final double MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;

    /**
     * Reads with length lower than this number, after clipping off overhangs outside the active region,
     * won't be considered for genotyping.
     */
    private static final int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    public static final int DEFAULT_READ_QUALITY_FILTER_THRESHOLD = 20;

    protected static final List<VariantContext> NO_CALLS = Collections.emptyList();

    private static final Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private static final Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    /**
     * Create and initialize a new HaplotypeCallerEngine given a collection of HaplotypeCaller arguments, a reads header,
     * and a reference file
     * @param hcArgs command-line arguments for the HaplotypeCaller
     * @param assemblyRegionArgs
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param readsHeader header for the reads
     * @param referenceReader reader to provide reference data, this reference reader will be closed when {@link #shutdown()} is called
     * @param annotationEngine variantAnnotatorEngine with annotations to process already added
     */
    public HaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, AssemblyRegionArgumentCollection assemblyRegionArgs, boolean createBamOutIndex,
                                 boolean createBamOutMD5, final SAMFileHeader readsHeader,
                                 CachingIndexedFastaSequenceFile referenceReader, VariantAnnotatorEngine annotationEngine) {
        this.dragstrParams = DragstrParamUtils.parse(hcArgs.likelihoodArgs.dragstrParams);
        this.hcArgs = Utils.nonNull(hcArgs);
        this.readsHeader = Utils.nonNull(readsHeader);
        this.referenceReader = Utils.nonNull(referenceReader);
        this.annotationEngine = Utils.nonNull(annotationEngine);
        this.aligner = SmithWatermanAligner.getAligner(hcArgs.smithWatermanImplementation);
        forceCallingAllelesPresent = hcArgs.alleles != null;

        // Add necessary debug streams to the output
        if (hcArgs.assemblyStateOutput != null) {
            try {
                assemblyDebugOutStream = new OutputStreamWriter(hcArgs.assemblyStateOutput.getOutputStream());
            } catch (final Exception e) {
                throw new UserException.CouldNotCreateOutputFile(hcArgs.assemblyStateOutput, "Provided argument for assembly debug graph location could not be created", e);
            }
        } else {
            assemblyDebugOutStream = null;
        }
        if (hcArgs.genotyperDebugOutStream != null) {
            HaplotypeCallerGenotypingDebugger.initialize(hcArgs.genotyperDebugOutStream);
        }

        // Parse the user provided custom ploidy regions into ploidyRegions object containing SimpleCounts
        if (this.hcArgs.ploidyRegions != null) {
            FeatureDataSource<NamedFeature> ploidyDataSource = new FeatureDataSource<>(this.hcArgs.ploidyRegions, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, NamedFeature.class);
            ploidyDataSource.forEach(r -> this.ploidyRegions.add(new SimpleCount(r)));
        }

        for (SimpleCount region : this.ploidyRegions) {
            if (!IntervalUtils.intervalIsOnDictionaryContig(region.getInterval(), readsHeader.getSequenceDictionary())) {
                throw new UserException("Invalid region provided for --ploidy-regions at " + region.getContig() + ":" + region.getStart() + "-" + region.getEnd() + ". Contig name or endpoint doesn't match read sequence dictionary.");
            }
        }

        this.ploidyRegionsOverlapDetector = OverlapDetector.create(this.ploidyRegions);
        this.allCustomPloidies = this.ploidyRegions.stream().map(SimpleCount::getCount).collect(Collectors.toCollection(LinkedHashSet::new));

        trimmer = new AssemblyRegionTrimmer(assemblyRegionArgs, readsHeader.getSequenceDictionary());
        initialize(createBamOutIndex, createBamOutMD5);
    }

    /**
     * Common method to use in order to remove unwanted annotations from the list returned by the plugin specifically
     * for reference confidence mode. Will also ensure StrandBiasBySample is present regardless of user requests.
     *
     * @param annotations a list of annotations to change
     * @return a list of annotations with non GVCF annotations removed
     */
    public static Collection<Annotation> filterReferenceConfidenceAnnotations(Collection<Annotation> annotations) {
        logger.info("Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled");

        // Override user preferences and add StrandBiasBySample
        if (!annotations.contains(new StrandBiasBySample())) {
            annotations.add(new StrandBiasBySample());
        }
        // Override user preferences and remove ChromosomeCounts, FisherStrand, StrandOddsRatio, and QualByDepth Annotations
        return annotations.stream()
                .filter(c -> (
                        c.getClass() != (ChromosomeCounts.class) &&
                        c.getClass() != (FisherStrand.class) &&
                        c.getClass() != (StrandOddsRatio.class) &&
                        c.getClass() != (QualByDepth.class))
                ).collect(Collectors.toList());
    }

    /**
     * @return the default set of variant annotations for use with HaplotypeCaller
     */
    public static List<Class<? extends Annotation>> getStandardHaplotypeCallerAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, StandardHCAnnotation.class);
    }

    private void initialize(boolean createBamOutIndex, final boolean createBamOutMD5) {
        // Note: order of operations matters here!

        initializeSamples();

        // Must be called after initializeSamples()
        validateAndInitializeArgs();
        minTailQuality = (byte)(hcArgs.minBaseQualityScore - 1);

        initializeActiveRegionEvaluationGenotyperEngine();

        defaultGenotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, ! hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
        defaultGenotypingEngine.setAnnotationEngine(annotationEngine);
        minosInitGenotypingVariants();
        minosInitWarm();
        minosInitRepair();

        // Create other custom genotyping engines if user provided ploidyRegions
        for (final int ploidy : this.allCustomPloidies) {
            HaplotypeCallerArgumentCollection newPloidyHcArgs = hcArgs.copyWithNewPloidy(ploidy);
            HaplotypeCallerGenotypingEngine newGenotypingEngine = new HaplotypeCallerGenotypingEngine(newPloidyHcArgs, samplesList, ! hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
            newGenotypingEngine.setAnnotationEngine(annotationEngine);
            this.ploidyToGenotyperMap.put(ploidy, newGenotypingEngine);
        }

        boolean isFlowBased = (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBased)
                || (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBasedHMM);
        referenceConfidenceModel = new ReferenceConfidenceModel(samplesList, readsHeader,
                hcArgs.indelSizeToEliminateInRefModel,
                hcArgs.standardArgs.genotypeArgs.numRefIfMissing,
                hcArgs.refModelDelQual,
                !hcArgs.overrideSoftclipFragmentCheck,
                isFlowBased);

        //Allele-specific annotations are not yet supported in the VCF mode
        if (isAlleleSpecificExceptHmerLengthOrStrandBiasMode(annotationEngine) && isVCFMode()){
           throw new UserException("Allele-specific annotations are not yet supported in the VCF mode");
        }

        haplotypeBAMWriter = AssemblyBasedCallerUtils.createBamWriter(hcArgs, createBamOutIndex, createBamOutMD5, readsHeader);
        alleleLikelihoodWriter = AssemblyBasedCallerUtils.createAlleleLikelihoodWriter(hcArgs);
        assemblyEngine = hcArgs.createReadThreadingAssembler();
        assembledEventMapVcfOutputWriter = Optional.ofNullable(hcArgs.assemblerArgs.debugAssemblyVariantsOut != null ?
                GATKVariantContextUtils.createVCFWriter(
                        new GATKPath(hcArgs.assemblerArgs.debugAssemblyVariantsOut).toPath(),
                        readsHeader.getSequenceDictionary(),
                        false,
                        Options.DO_NOT_WRITE_GENOTYPES, Options.INDEX_ON_THE_FLY)
                : null);
        assembledEventMapVariants = Optional.ofNullable(hcArgs.assemblerArgs.debugAssemblyVariantsOut != null ?
                new PriorityQueue<>(200, new VariantContextComparator(readsHeader.getSequenceDictionary())) : null);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> writeHeader(writer, readsHeader.getSequenceDictionary(), new HashSet<>()));
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs, !hcArgs.softClipLowQualityEnds, hcArgs.likelihoodArgs.likelihoodEngineImplementation);


        //Some sanity checking of the stepwise filtering approach
        if (hcArgs.stepwiseFiltering) {
            if (hcArgs.likelihoodArgs.likelihoodEngineImplementation == ReadLikelihoodCalculationEngine.Implementation.FlowBased) {
                throw new UserException("'--"+ HaplotypeCallerArgumentCollection.STEPWISE_FITLERING_ARGUMENT + "' is specified but the selected likelihoods engine is FlowBased. This is invalid.");
            }
            if (!hcArgs.filterAlleles) {
                throw new UserException("'--"+ HaplotypeCallerArgumentCollection.STEPWISE_FITLERING_ARGUMENT + "' is specified, but allele filtering is disabled. Please enable allele filtering with '--"+AssemblyBasedCallerArgumentCollection.FILTER_ALLELES+"'");
            }
        }

        filterStepLikelihoodCalculationEngine = (hcArgs.stepwiseFiltering?
                AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs, !hcArgs.softClipLowQualityEnds, ReadLikelihoodCalculationEngine.Implementation.FlowBased)
                : null);
        pdhmmLikelihoodCalculationEngine = (hcArgs.pileupDetectionArgs.generatePDHaplotypes?
                new PDPairHMMLikelihoodCalculationEngine((byte) hcArgs.likelihoodArgs.gcpHMM, hcArgs.likelihoodArgs.dontUseDragstrPairHMMScores ? null : DragstrParamUtils.parse(hcArgs.likelihoodArgs.dragstrParams),
                        hcArgs.likelihoodArgs.pairHMMNativeArgs.getPairHMMArgs(), hcArgs.likelihoodArgs.pairHMM, hcArgs.pileupDetectionArgs.pdhmmDebugOutputResults, AssemblyBasedCallerUtils.getGlobalMismatchingRateFromArgs(hcArgs.likelihoodArgs), hcArgs.likelihoodArgs.pcrErrorModel,
                        hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD, hcArgs.likelihoodArgs.enableDynamicReadDisqualification, hcArgs.likelihoodArgs.readDisqualificationThresholdConstant,
                        hcArgs.likelihoodArgs.expectedErrorRatePerBase, !hcArgs.likelihoodArgs.disableSymmetricallyNormalizeAllelesToReference, hcArgs.likelihoodArgs.disableCapReadQualitiesToMapQ, !hcArgs.softClipLowQualityEnds,
                        hcArgs.pileupDetectionArgs.pdhmmOptimization? hcArgs.informativeReadOverlapMargin : -1) //Logic to control whether we apply the pdhmm optimization
                : null);
    }

    private boolean isVCFMode() {
        return hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.NONE;
    }

    private boolean isAlleleSpecificExceptHmerLengthOrStrandBiasMode(final VariantAnnotatorEngine annotationEngine) {
        //HACK. Note: we can't use subclass information from ReducibleAnnotation (which would be the obvious choice)
        // because RMSMappingQuality is both a reducible annotation and a standard annotation.

        return annotationEngine.getInfoAnnotations().stream()
                .filter(infoFieldAnnotation -> !infoFieldAnnotation.getClass().getSimpleName().equals(AS_StrandBiasMutectAnnotation.class.getSimpleName()))
                .anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().startsWith("AS_")) ||
                annotationEngine.getGenotypeAnnotations().stream()
                .anyMatch(genotypeAnnotation -> genotypeAnnotation.getClass().getSimpleName().startsWith("AS_"));
    }

    private void validateAndInitializeArgs() {
        if ( hcArgs.standardArgs.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is only available for diploid samples.");
        }

        if ( hcArgs.dontGenotype && emitReferenceConfidence() ) {
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");
        }

        if ( emitReferenceConfidence() ) {
            hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling = -0.0;
            logger.info("Standard Emitting and Calling confidence set to " + hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling + " for reference-model confidence output");
            if ( ! hcArgs.standardArgs.annotateAllSitesWithPLs ) {
                logger.info("All sites annotated with PLs forced to true for reference-model confidence output");
            }
            hcArgs.standardArgs.annotateAllSitesWithPLs = true;
        } else if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            hcArgs.doNotRunPhysicalPhasing = true;
            logger.info("Disabling physical phasing, which is supported only for reference-model confidence output");
        }

        if ( hcArgs.standardArgs.CONTAMINATION_FRACTION_FILE != null ) {
            hcArgs.standardArgs.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(hcArgs.standardArgs.CONTAMINATION_FRACTION_FILE, hcArgs.standardArgs.CONTAMINATION_FRACTION, sampleSet, logger));
        }

        Utils.validateArg(hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD >= QualityUtils.MIN_USABLE_Q_SCORE, "BASE_QUALITY_SCORE_THRESHOLD must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");

        if ( emitReferenceConfidence() && samplesList.numberOfSamples() != 1 ) {
            throw new CommandLineException.BadArgumentValue(AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, "Can only be used in single sample mode currently. Use the --sample-name argument to run on a single sample out of a multi-sample BAM file.");
        }

        if (hcArgs.floorBlocks && !emitReferenceConfidence()) {
            throw new UserException(HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS + " refers to GVCF blocks," +
                    " so reference confidence mode (" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME +
                    ") must be specified.");
        }
    }

    private void initializeSamples() {
        sampleSet = ReadUtils.getSamplesFromHeader(readsHeader);
        samplesList = new IndexedSampleList(sampleSet);

        if ( hcArgs.sampleNameToUse != null ) {
            if ( ! sampleSet.contains(hcArgs.sampleNameToUse) ) {
                throw new CommandLineException.BadArgumentValue("--sample_name", "Specified name does not exist in input bam files");
            }
            if (sampleSet.size() == 1) {
                //No reason to incur performance penalty associated with filtering if they specified the name of the only sample
                hcArgs.sampleNameToUse = null;
            } else {
                samplesList = new IndexedSampleList(hcArgs.sampleNameToUse);
                sampleSet.clear();
                sampleSet.add(hcArgs.sampleNameToUse);
            }
        }
    }

    private void initializeActiveRegionEvaluationGenotyperEngine() {
        final StandardCallerArgumentCollection activeRegionArgs = new StandardCallerArgumentCollection();
        activeRegionArgs.copyStandardCallerArgsFrom(hcArgs.standardArgs);

        activeRegionArgs.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        activeRegionArgs.genotypeArgs.standardConfidenceForCalling = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, hcArgs.standardArgs.genotypeArgs.standardConfidenceForCalling); // low values used for isActive determination only, default/user-specified values used for actual calling
        activeRegionArgs.CONTAMINATION_FRACTION = 0.0;
        activeRegionArgs.CONTAMINATION_FRACTION_FILE = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use ploidy == 1
        activeRegionArgs.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, hcArgs.standardArgs.genotypeArgs.samplePloidy);

        defaultActiveRegionEvaluationGenotyperEngine = new MinimalGenotypingEngine(activeRegionArgs, samplesList);
        defaultActiveRegionEvaluationGenotyperEngine.setLogger(logger);

        // If custom ploidyRegions provided, create corresponding map for active region determination genotyper
        for (final int ploidy : this.allCustomPloidies) {
            StandardCallerArgumentCollection newPloidyActiveArgs = new StandardCallerArgumentCollection();
            newPloidyActiveArgs.copyStandardCallerArgsFrom(activeRegionArgs);
            newPloidyActiveArgs.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, ploidy);
            MinimalGenotypingEngine newActiveGenotypingEngine = new MinimalGenotypingEngine(newPloidyActiveArgs, samplesList);
            newActiveGenotypingEngine.setLogger(logger);
            this.ploidyToActiveEvaluationGenotyper.put(ploidy, newActiveGenotypingEngine);
        }
    }

    /**
     * @return the default set of read filters for use with the HaplotypeCaller
     */
    public static List<ReadFilter> makeStandardHCReadFilters() {
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(DEFAULT_READ_QUALITY_FILTER_THRESHOLD));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }

    /**
     * @return the default read transformer to apply to HaplotypeCaller input
     */
    public static ReadTransformer makeStandardHCReadTransformer() {
        return new IUPACReadTransformer(true);
    }

    /**
     * Create a VCF or GVCF writer as appropriate, given our arguments
     *
     * @param outputVCF location to which the vcf should be written
     * @param readsDictionary sequence dictionary for the reads
     * @return a VCF or GVCF writer as appropriate, ready to use
     */
    public VariantContextWriter makeVCFWriter( final GATKPath outputVCF, final SAMSequenceDictionary readsDictionary,
                                               final boolean createOutputVariantIndex, final boolean  createOutputVariantMD5,
                                               final boolean sitesOnlyMode ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        final List<Options> options = new ArrayList<>(2);
        if (createOutputVariantIndex) {options.add(Options.INDEX_ON_THE_FLY);}
        if (sitesOnlyMode) {options.add(Options.DO_NOT_WRITE_GENOTYPES);}

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputVCF.toPath(),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

        if ( hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
            try {
                writer = new GVCFWriter(writer, hcArgs.GVCFGQBands, hcArgs.floorBlocks);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    /**
     * Create a VCF header.
     *
     * @param sequenceDictionary sequence dictionary for the reads
     * @return a VCF header
     */
    public VCFHeader makeVCFHeader( final SAMSequenceDictionary sequenceDictionary, final Set<VCFHeaderLine>  defaultToolHeaderLines ) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(defaultToolHeaderLines);

        headerInfo.addAll(defaultGenotypingEngine.getAppropriateVCFInfoHeaders());
        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions(emitReferenceConfidence()));
        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        if (hcArgs.standardArgs.genotypeArgs.supportVariants != null) {
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
            headerInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        }

        if ( ! hcArgs.doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
            headerInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.PHASE_SET_KEY));
        }

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        if (( hcArgs.flowAssemblyCollapseHKerSize > 0 ) || (hcArgs.flowAssemblyCollapseHKerSize==AssemblyBasedCallerUtils.DETERMINE_COLLAPSE_THRESHOLD)){
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXT_COLLAPSED_KEY));
        }

        if (hcArgs.filterAlleles) {
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSSIBLE_FP_ADJACENT_TP_KEY));
        }

        if ( emitReferenceConfidence() ) {
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
        }

        if (hcArgs.likelihoodArgs.dragstrParams != null) {
            headerInfo.addAll(DragstrVariantContextAnnotations.vcfHeaderLines());
        }

        if (hcArgs.standardArgs.genotypeArgs.genotypeAssignmentMethod == GenotypeAssignmentMethod.USE_POSTERIOR_PROBABILITIES) {
            headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Float, "genotype posterior in Phred Scale" ));
            headerInfo.add(new VCFFormatHeaderLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Float, "genotype priors in Phred Scale"));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo, sampleSet);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        return vcfHeader;
    }

    private SAMSequenceDictionary minosSequenceDictionary;
    private Set<VCFHeaderLine> minosDefaultHeaderLines;

    /**
     * Writes an appropriate VCF header, given our arguments, to the provided writer
     *
     * @param vcfWriter writer to which the header should be written
     */
    public void writeHeader( final VariantContextWriter vcfWriter, final SAMSequenceDictionary sequenceDictionary,
                             final Set<VCFHeaderLine>  defaultToolHeaderLines) {
        Utils.nonNull(vcfWriter);
        minosSequenceDictionary = sequenceDictionary;
        minosDefaultHeaderLines = defaultToolHeaderLines;
        vcfWriter.writeHeader(makeVCFHeader(sequenceDictionary, defaultToolHeaderLines));
        for (final MinosGenotypingVariant v : minosGenotypingVariants) {
            v.writer = GATKVariantContextUtils.createVCFWriter(Paths.get(v.output), sequenceDictionary, false,
                    htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY);
            v.writer.writeHeader(makeVCFHeader(sequenceDictionary, defaultToolHeaderLines));
        }
    }

    /**
     * Determines the appropriate ploidy to use at given site for different genotyping engines.
     * @param region Current region of interest
     * @return Ploidy value to use here given user inputs, or -1 if fall back to default
     */
    private int getPloidyToUseAtThisSite(Locatable region) {
        Set<SimpleCount> overlaps = this.ploidyRegionsOverlapDetector.getOverlaps(region);
        // Return first engine for interval overlapping this region
        if (!overlaps.isEmpty()) {
            Iterator<SimpleCount> intervals = overlaps.iterator();
            int max = intervals.next().getCount();
            while (intervals.hasNext()) {
                int next = intervals.next().getCount();
                if (next > max) {
                    max = next;
                }
            }
            return max;
        } else {
            return -1; // Sentinel value to fall back to default genotyper
        }
    }

    /**
     * Selects appropriate active region genotyping engine for given region
     * @param region Current region of interest
     * @return Active genotyping engine with correct ploidy setting for given region
     */
    private MinimalGenotypingEngine getLocalActiveGenotyper(Locatable region) {
        int currentPloidy = getPloidyToUseAtThisSite(region);
        if (currentPloidy == -1) {
            return this.defaultActiveRegionEvaluationGenotyperEngine;
        } else {
            return this.ploidyToActiveEvaluationGenotyper.get(currentPloidy);
        }
    }

    /**
     * Selects appropriate genotyping engine for given region.
     * @param region Current region of interest, e.g. AssemblyRegion
     * @return Genotyping engine with correct ploidy setting for given region
     */
    protected HaplotypeCallerGenotypingEngine getLocalGenotypingEngine(Locatable region) {
        int currentPloidy = getPloidyToUseAtThisSite(region);
        if (currentPloidy == -1) {
            return this.defaultGenotypingEngine;
        } else {
            return this.ploidyToGenotyperMap.get(currentPloidy);
        }
    }

    /**
     * Given a pileup, returns an ActivityProfileState containing the probability (0.0 to 1.0) that it's an "active" site.
     *
     * Note that the current implementation will always return either 1.0 or 0.0, as it relies on the smoothing in
     * {@link org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile} to create the full distribution
     * of probabilities. This is consistent with GATK 3.
     *
     * @param context reads pileup to examine
     * @param ref reference base overlapping the pileup locus
     * @param features features overlapping the pileup locus
     * @return probability between 0.0 and 1.0 that the site is active (in practice with this implementation: either 0.0 or 1.0)
     */
    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext features) {
        MinimalGenotypingEngine localActiveGenotypingEngine = getLocalActiveGenotyper(ref);

        if (forceCallingAllelesPresent && features.getValues(hcArgs.alleles, ref).stream().anyMatch(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())) {
            return new ActivityProfileState(ref.getInterval(), 1.0);
        }

        if (context == null || context.getBasePileup().isEmpty()) {
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getInterval(), 0.0);
        }

        final int ploidy = localActiveGenotypingEngine.getConfiguration().genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied

        final Map<String, AlignmentContext> splitContexts;
        if (samplesList.numberOfSamples() == 1) {
            //If we know a priori that there's just one sample, take a shortcut and dont examine each read in the pileup
            splitContexts = context.splitContextBySampleName(samplesList.getSample(0), readsHeader);
        } else {
            splitContexts = context.splitContextBySampleName(readsHeader);
        }

        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        for (final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet()) {
            // We summarize the reads mismatches vs the ref at this stage to avoid issues due to reads being hardclipped down the road.
            if (hcArgs.pileupDetectionArgs.usePileupDetection) {
                sample.getValue().getBasePileup().forEach(p -> PileupBasedAlleles.addMismatchPercentageToRead(p.getRead(), readsHeader, ref));
            }
            // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
            final int activeRegionDetectionHackishSamplePloidy = localActiveGenotypingEngine.getConfiguration().genotypeArgs.samplePloidy;
            final double[] genotypeLikelihoods = ((RefVsAnyResult) referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(
                    activeRegionDetectionHackishSamplePloidy,
                    sample.getValue().getBasePileup(), ref.getBase(),
                    hcArgs.minBaseQualityScore,
                    averageHQSoftClips, false)).genotypeLikelihoods;
            genotypes.add(new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make());
        }

        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE, FAKE_ALT_ALLELE);
        final double isActiveProb;

        if (genotypes.size() == 1) {
            // Faster implementation using exact marginalization instead of iteration
            isActiveProb = localActiveGenotypingEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());
        } else {
            final VariantContext vcOut = localActiveGenotypingEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getEnd(), alleles).genotypes(genotypes).make());
            isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());
        }

        // Here we store the crude fast-likelihood score used to determine positional activity for the ActivityProfileState.
        // We need to use this score later for DRAGEN-GATK, as one of the Heuristics applied by the new Pileup-Caller is that
        // now we optionally will ignore pileup events from positions on the genome that did not ping as being "active," thus
        // we must take care to preserve this information for later.
        //TODO If you dare, the entire ##PileupBasedAlleles.getPileupVariantContexts() step could be moved here to cut down on storing these intermediate scores for later
        //TODO This equivalent operation does not currently exist in the #Mutect2.isActive() method
        ActivityProfileState output = new ActivityProfileState(ref.getInterval(), isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean());
        double[] gts = genotypes.get(0).getLikelihoods().getAsVector();
        int maxElementIndex = MathUtils.maxElementIndex(gts);
        output.setOriginalActiveProb(gts[maxElementIndex] - gts[0]);
        return output;

    }

    /**
     * Generate variant calls for an assembly region
     *
     * @param region region to assemble and perform variant calling on
     * @param features Features overlapping the assembly region
     * @return List of variants discovered in the region (may be empty)
     */
    public List<VariantContext> callRegion(final AssemblyRegion region, final FeatureContext features, final ReferenceContext referenceContext) {
        if (minosRepairLog != null) {
            return minosRepairRegion(region, features, referenceContext);
        }
        if (minosWarmDir != null) {
            minosWarmKey = minosRegionKey(region);
        }
        if (minosGenotypingVariants.isEmpty()) {
            return minosCallRegion(region, features, referenceContext);
        }
        // minos: log what every main-region variant did in this region (read fingerprint + random draws), so a
        // diverged variant can later be repaired region by region
        final String key = minosRegionKey(region);
        final List<MinosGenotypingVariant> logged = new ArrayList<>();
        final List<long[]> before = new ArrayList<>();
        for (final MinosGenotypingVariant v : minosGenotypingVariants) {
            if (v.stream < 0) {
                logged.add(v);
                before.add(new long[]{v.rng.gaussians(), v.well.steps()});
            }
        }
        try {
            return minosCallRegion(region, features, referenceContext);
        } finally {
            for (int i = 0; i < logged.size(); i++) {
                logged.get(i).logRegion(key, before.get(i));
            }
        }
    }

    private List<VariantContext> minosCallRegion(final AssemblyRegion region, final FeatureContext features, final ReferenceContext referenceContext) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(region);

        if ( hcArgs.justDetermineActiveRegions ) {
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;
        }
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("calling for region: " +region.getSpan());
        }

        final List<VariantContext> VCpriors = new ArrayList<>();
        if (hcArgs.standardArgs.genotypeArgs.supportVariants != null) {
            features.getValues(hcArgs.standardArgs.genotypeArgs.supportVariants).stream().forEach(VCpriors::add);
        }

        if ( hcArgs.sampleNameToUse != null ) {
            removeReadsFromAllSamplesExcept(hcArgs.sampleNameToUse, region);
        }

        if( ! region.isActive() ) {
            // Not active so nothing to do!
            return referenceModelForNoVariation(region, true, VCpriors);
        }

        final List<Event> givenAlleles = features.getValues(hcArgs.alleles).stream()
                .filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())
                .flatMap(vc -> GATKVariantContextUtils.splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
                .collect(Collectors.toList());

        if( givenAlleles.isEmpty() && region.size() == 0 ) {
            // No reads here so nothing to do!
            return referenceModelForNoVariation(region, true, VCpriors);
        }

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\n\n\n\n" + region.getSpan() + "\nNumber of reads in region: " + region.getReads().size() + "     they are:\n");
                for (GATKRead read : region.getReads()) {
                    assemblyDebugOutStream.write(read.getName() + "   " + read.convertToSAMRecord(region.getHeader()).getFlags() + "\n");
                }
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(region, hcArgs, readsHeader, samplesList, logger, referenceReader, assemblyEngine, aligner, !hcArgs.doNotCorrectOverlappingBaseQualities, hcArgs.fbargs, false);
        minosCallAssemblyVariants(region, features, referenceContext, givenAlleles);
        ReadThreadingAssembler.addAssembledVariantsToEventMapOutput(untrimmedAssemblyResult, assembledEventMapVariants, hcArgs.maxMnpDistance, assembledEventMapVcfOutputWriter);

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\nThere were " + untrimmedAssemblyResult.getHaplotypeList().size() + " haplotypes found. Here they are:\n");
                for (final String haplotype : untrimmedAssemblyResult.getHaplotypeList().stream().map(Haplotype::toString).sorted().collect(Collectors.toList())) {
                    assemblyDebugOutStream.write(haplotype);
                    assemblyDebugOutStream.append('\n');
                }
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        Pair<Set<Event>, Set<Event>> goodAndBadPileupEvents =
                PileupBasedAlleles.goodAndBadPileupEvents(region.getAlignmentData(), hcArgs.pileupDetectionArgs, readsHeader, hcArgs.minBaseQualityScore);
        final Set<Event> goodPileupEvents = goodAndBadPileupEvents.getLeft();
        final Set<Event> badPileupEvents = goodAndBadPileupEvents.getRight();

        // Regenerate the list of AllVariationEvents, filtering out assembled variants that must be filtered according to the pileupcaller code.
        final SortedSet<Event> allVariationEvents = untrimmedAssemblyResult.getVariationEvents(hcArgs.maxMnpDistance).stream()
                .filter(event -> !badPileupEvents.contains(event))
                .collect(Collectors.toCollection(() -> new TreeSet<>(AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR)));

        goodPileupEvents.forEach(allVariationEvents::add);
        givenAlleles.forEach(allVariationEvents::add);

        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(region, allVariationEvents, referenceContext);

        if (!trimmingResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        AssemblyResultSet assemblyResult = untrimmedAssemblyResult.trimTo(trimmingResult.getVariantRegion());
        assemblyResult.addGivenAlleles(givenAlleles, hcArgs.maxMnpDistance, aligner, hcArgs.getHaplotypeToReferenceSWParameters());

        // Pre-work for the PDHMM, if we are in PDHMM mode then here is where we re-compute the haplotypes as PD haplotypes.
        if (hcArgs.pileupDetectionArgs.generatePDHaplotypes) {
            // Note: we ignore maxMNPDistance here because the PartiallyDeterminedHaplotypeComputationEngine does not currently handle MNPs
            SortedSet<Event> assemblyVariants = assemblyResult.getVariationEvents(0);
            if (hcArgs.pileupDetectionArgs.debugPileupStdout) {
                System.out.println("Generating PDHaplotypes for PDHMM");
                System.out.println("Assembled Variants to use:");
                assemblyVariants.forEach(System.out::println);
                System.out.println("Good Pileup Variants to use:");
                goodPileupEvents.forEach(System.out::println);
                System.out.println("Bad Pileup Variants to filter:");
                badPileupEvents.forEach(System.out::println);
                System.out.println("Adding Variants To Reference Haplotype:");
                System.out.println(assemblyResult.getReferenceHaplotype());
                System.out.println("FinalSpan: " + assemblyResult.getReferenceHaplotype().getGenomeLocation());
                System.out.println("CallingSpan: " + region.getSpan());
            }
            assemblyResult = PartiallyDeterminedHaplotypeComputationEngine.generatePDHaplotypes(assemblyResult,
                    badPileupEvents,
                    goodPileupEvents,
                    aligner,
                    hcArgs);
        }

        // Legacy Pileupcaller code. Supplement the assembly haps with artifical haps constructed from the discovered pileupcaller
        // alleles based off of the GenotypeGivenAlleles code approach.
        // NOTE: We might also call this if hcArgs.pileupDetectionArgs.useGGAFallback is set and we are making PD haplotypes. This
        //       fallback handles edge cases where the PartiallyDeterminedHaplotypeComputationEngine generates too-many haplotypes
        //       from a very complex assebmly region and we want to fall back to assembly.
        if (!hcArgs.pileupDetectionArgs.generatePDHaplotypes ||
                        (hcArgs.pileupDetectionArgs.useGGAFallback && !assemblyResult.hasOverwrittenHaps())) { // If we are generating PDHaps assert that it failed before calling this
            assemblyResult.removeHaplotypesWithBadAlleles(hcArgs, badPileupEvents);
            assemblyResult.injectPileupEvents(region, hcArgs, aligner, goodPileupEvents);
        }
        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        final List<GATKRead> readStubs = regionForGenotyping.getReads().stream()
                .filter(r -> AlignmentUtils.unclippedReadLength(r)  < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        regionForGenotyping.removeAll(readStubs);

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKRead>> perSampleFilteredReadList = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, filteredReads);

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (!assemblyResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        // For sure this is not true if gVCF is on.
        if (hcArgs.dontGenotype) {
            return NO_CALLS; // user requested we not proceed
        }

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (regionForGenotyping.size() == 0 && !hcArgs.disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(region, false, VCpriors);
        }

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("\n=================================================");
            HaplotypeCallerGenotypingDebugger.println("assemblyRegion: "+new SimpleInterval(region));
            HaplotypeCallerGenotypingDebugger.println("=================================================");
        }

        // evaluate each sample's reads against all haplotypes
        final Map<String,List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, regionForGenotyping.getReads());
        final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods;
        List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();


        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("\nUnclipped Haplotypes("+haplotypes.size()+"):");
            for (Haplotype haplotype : untrimmedAssemblyResult.getHaplotypeList()) {
                HaplotypeCallerGenotypingDebugger.println("["+haplotype.getStart()+"-"+haplotype.getEnd()+"] k="+haplotype.getKmerSize()+" len: "+haplotype.length()+" "+haplotype.getCigar()+(haplotype.isReference()?"ref":""));
                HaplotypeCallerGenotypingDebugger.println(haplotype.toString());
            }

            HaplotypeCallerGenotypingDebugger.println("\nClipped Haplotyes("+haplotypes.size()+"):");
            for (Haplotype haplotype : haplotypes) {
                HaplotypeCallerGenotypingDebugger.println("["+haplotype.getStart()+"-"+haplotype.getEnd()+"] k="+haplotype.getKmerSize()+" len: "+haplotype.length()+" "+haplotype.getCigar()+(haplotype.isReference()?"ref":""));
                HaplotypeCallerGenotypingDebugger.println(haplotype.toString());
            }
            HaplotypeCallerGenotypingDebugger.println("");
        }

        // Calculate the likelihoods: CPU intensive part.
        // flow based alignment might add an extra step of uncollapsing - implemented by possiblyUncollapseHaplotypesInReadLikelihoods
        // non-flow based alignment will not be affected.
        readLikelihoods = possiblyUncollapseHaplotypesInReadLikelihoods(untrimmedAssemblyResult,
                hcArgs.stepwiseFiltering
                        ? filterStepLikelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true)
                        : ( assemblyResult.isPartiallyDeterminedList() ?
                        pdhmmLikelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true) :
                        likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads, true)));

        alleleLikelihoodWriter.ifPresent(
                writer -> writer.writeAlleleLikelihoods(readLikelihoods));

        // Optional: pre-filter haplotypes by removing weak/noisy alleles
        // this is important for the genotyping, as weak allele close to the real allele often decreases its quality
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;
        Set<Integer> suspiciousLocations = new HashSet<>();
        if (hcArgs.filterAlleles) {
            logger.debug("Filtering alleles");

            AlleleFilteringHC alleleFilter = new AlleleFilteringHC(hcArgs, assemblyDebugOutStream, localGenotypingEngine);
            //need to update haplotypes to find the alleles
            EventMap.buildEventMapsForHaplotypes(readLikelihoods.alleles(),
                    assemblyResult.getFullReferenceWithPadding(),
                    assemblyResult.getPaddedReferenceLoc(),
                    hcArgs.assemblerArgs.debugAssembly,
                    hcArgs.maxMnpDistance);
            subsettedReadLikelihoodsFinal = alleleFilter.filterAlleles(readLikelihoods, assemblyResult.getPaddedReferenceLoc().getStart(), suspiciousLocations);

        } else {
            subsettedReadLikelihoodsFinal = readLikelihoods;
        }


        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]
        haplotypes = subsettedReadLikelihoodsFinal.alleles();

        // Reproscess with the HMM if we are in stepwise filtering mode
        if (hcArgs.stepwiseFiltering) {
            subsettedReadLikelihoodsFinal = (likelihoodCalculationEngine).computeReadLikelihoods(
                    subsettedReadLikelihoodsFinal.alleles(),
                    assemblyResult.getRegionForGenotyping().getHeader(),
                    samplesList, reads, true);
        }

        //Realign reads to their best haplotype.
        final SWParameters readToHaplotypeSWParameters = hcArgs.getReadToHaplotypeSWParameters();
        // TODO Yes we skip realignment entirely when we are in DRAGEN-GATK PDHMM mode. Realignment of the reads makes no sense when
        // TODO the bases of the haplotypes used for calling no longer reflect specified variants present.
        if (!(hcArgs.pileupDetectionArgs.generatePDHaplotypes && !hcArgs.pileupDetectionArgs.useDeterminedHaplotypesDespitePdhmmMode)) {
            final Map<GATKRead, GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(subsettedReadLikelihoodsFinal, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner, readToHaplotypeSWParameters);
            subsettedReadLikelihoodsFinal.changeEvidence(readRealignments);
        }

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            for (int counter = 0; counter < readLikelihoods.sampleEvidence(0).size(); counter++) {
                GATKRead read = readLikelihoods.sampleEvidence(0).get(counter);
                HaplotypeCallerGenotypingDebugger.println("read " + counter + ": " + read.getName() + " cigar: " + read.getCigar() + " mapQ: " + read.getMappingQuality() + " loc: [" + read.getStart() + "-" + read.getEnd() + "] unclippedloc: [" + read.getUnclippedStart() + "-" + read.getUnclippedEnd() + "] length:" + read.getLength());
                String hmmScores = "";
                for (int h = 0; h < readLikelihoods.numberOfAlleles(); h++) {
                    hmmScores = hmmScores + "," + readLikelihoods.sampleMatrix(0).get(h, counter);
                }
                HaplotypeCallerGenotypingDebugger.println(hmmScores);
            }
        }

        minosBeforeGenotyping();
        final CalledHaplotypes calledHaplotypes = localGenotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                subsettedReadLikelihoodsFinal,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getSpan(),
                features,
                givenAlleles,
                emitReferenceConfidence(),
                hcArgs.maxMnpDistance,
                readsHeader,
                haplotypeBAMWriter.isPresent(),
                suspiciousLocations,
                readLikelihoods);
        // minos: one extra genotyping pass per variant. Variants whose PairHMM parameters equal the main config's
        // reuse its likelihoods; the others share one likelihood computation per distinct PairHMM parameter set,
        // on the same region, haplotypes and reads (computeReadLikelihoods copies its inputs; the only thing it
        // writes on shared reads is a transient attribute read back within the same call). Each variant has its
        // own genotyping engine and random streams.
        final Map<String, AlleleLikelihoods<GATKRead, Haplotype>> minosLikelihoods = new HashMap<>();
        for (final MinosGenotypingVariant v : minosGenotypingVariants) {
            if (v.assemblyKey != null || v.stream >= 0) {
                continue;                          // handled right after the main assembly / by its activity stream
            }
            AlleleLikelihoods<GATKRead, Haplotype> rl = subsettedReadLikelihoodsFinal;
            AlleleLikelihoods<GATKRead, Haplotype> preFilter = readLikelihoods;
            if (v.likelihoodKey != null) {
                rl = minosLikelihoods.get(v.likelihoodKey);
                if (rl == null) {
                    rl = possiblyUncollapseHaplotypesInReadLikelihoods(untrimmedAssemblyResult,
                            minosLikelihoodEngines.get(v.likelihoodKey).computeReadLikelihoods(assemblyResult, samplesList, reads, true));
                    // exactly the main path's next step: genotyping sees reads realigned to their best haplotype
                    if (!(hcArgs.pileupDetectionArgs.generatePDHaplotypes && !hcArgs.pileupDetectionArgs.useDeterminedHaplotypesDespitePdhmmMode)) {
                        rl.changeEvidence(AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(rl,
                                assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner,
                                readToHaplotypeSWParameters));
                    }
                    minosLikelihoods.put(v.likelihoodKey, rl);
                }
                preFilter = rl;
            }
            Utils.minosSetRandomOverride(v.rng);          // its QD jitter comes from its own stream, as in its own run
            Utils.minosSetRandomDataOverride(v.rdg);      // and its contamination downsampling too
            try {
                v.engine.assignGenotypeLikelihoods(rl.alleles(), rl, perSampleFilteredReadList,
                        assemblyResult.getFullReferenceWithPadding(), assemblyResult.getPaddedReferenceLoc(),
                        regionForGenotyping.getSpan(), features, givenAlleles, emitReferenceConfidence(), hcArgs.maxMnpDistance,
                        readsHeader, haplotypeBAMWriter.isPresent(), suspiciousLocations, preFilter)
                    .getCalls().stream().map(RMSMappingQuality.getInstance()::finalizeRawMQ).forEach(v.writer::add);
            } finally {
                Utils.minosSetRandomOverride(null);
                Utils.minosSetRandomDataOverride(null);
            }
        }
        minosAfterGenotyping();
        if (minosWarmActive) {
            minosWarmCapture(subsettedReadLikelihoodsFinal, perSampleFilteredReadList,
                    assemblyResult.getFullReferenceWithPadding(), assemblyResult.getPaddedReferenceLoc(),
                    regionForGenotyping.getSpan(), features, givenAlleles, readsHeader, suspiciousLocations);
        }

        if ( haplotypeBAMWriter.isPresent() ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if ( hcArgs.disableOptimizations ) {
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            }
            haplotypeBAMWriter.get().writeReadsAlignedToHaplotypes(haplotypes, assemblyResult.getPaddedReferenceLoc(), haplotypes,
                                                             calledHaplotypeSet, subsettedReadLikelihoodsFinal,regionForGenotyping.getSpan());
        }

        if( hcArgs.assemblerArgs.debugAssembly) {
            logger.info("----------------------------------------------------------------------------------");
        }

        if ( emitReferenceConfidence() ) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(region, false, VCpriors);
            }
            else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section, then variant-containing section, then right flank
                trimmingResult.nonVariantLeftFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, VCpriors)));

                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        subsettedReadLikelihoodsFinal, localGenotypingEngine.getPloidyModel(), calledHaplotypes.getCalls(), hcArgs.standardArgs.genotypeArgs.supportVariants != null,
                        VCpriors));

                trimmingResult.nonVariantRightFlankRegion().ifPresent(flank -> result.addAll(referenceModelForNoVariation(flank, false, VCpriors)));

                return result;
            }
        }
        else {
            //TODO this should be updated once reducible annotations are handled properly.
            return calledHaplotypes.getCalls()
                    .stream()
                    .map(RMSMappingQuality.getInstance()::finalizeRawMQ)
                    .collect(Collectors.toList());
        }
    }

    /*
     * possibly uncollapse the haplotypes inside a read likelihood matrix
     *
     * At this stage, the Haplotypes in the likelihod matrix are derived from the reads. If the reads
     * are flow reads, then they are essentially collapsed, i.e. their maximal hmer size is limited
     * by the flow format's maxHmer. The method uncollapses the haplotypes to be consistant, as much
     * as possible (or applicable) according to the reference.
     */
    private AlleleLikelihoods<GATKRead, Haplotype> possiblyUncollapseHaplotypesInReadLikelihoods(final AssemblyResultSet assemblyResultSet, final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods) {

        if ( assemblyResultSet.getHaplotypeCollapsingEngine() != null ) {

            final List<Haplotype> haplotypes = assemblyResultSet.getHaplotypeCollapsingEngine().uncollapseHmersInHaplotypes(readLikelihoods.alleles(), false, null);
            logger.debug(() -> String.format("%d haplotypes before uncollapsing", haplotypes.size()));
            Map<Haplotype, List<Haplotype>> identicalHaplotypesMap = LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence(haplotypes);
            readLikelihoods.changeAlleles(haplotypes);
            final AlleleLikelihoods<GATKRead, Haplotype>  uncollapsedReadLikelihoods = readLikelihoods.marginalize(identicalHaplotypesMap);
            logger.debug(() -> String.format("%d haplotypes after uncollapsing", uncollapsedReadLikelihoods.numberOfAlleles()));

            return uncollapsedReadLikelihoods;
        } else {
            return readLikelihoods;
        }
    }


    protected boolean containsCalls(final CalledHaplotypes calledHaplotypes) {
        return calledHaplotypes.getCalls().stream()
                .flatMap(call -> call.getGenotypes().stream())
                .anyMatch(Genotype::isCalled);
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    protected List<VariantContext> referenceModelForNoVariation(final AssemblyRegion region, final boolean needsToBeFinalized, final List<VariantContext> VCpriors) {
        final HaplotypeCallerGenotypingEngine localGenotypingEngine = getLocalGenotypingEngine(region);
        if ( emitReferenceConfidence() ) {
            if ( needsToBeFinalized ) {
                AssemblyBasedCallerUtils.finalizeRegion(region,
                        hcArgs.assemblerArgs.errorCorrectReads,
                        hcArgs.dontUseSoftClippedBases,
                        minTailQuality, readsHeader, samplesList,
                        ! hcArgs.doNotCorrectOverlappingBaseQualities,
                        hcArgs.softClipLowQualityEnds,
                        hcArgs.overrideSoftclipFragmentCheck,
                        hcArgs.pileupDetectionArgs.usePileupDetection);
            }

            filterNonPassingReads(region);

            final SimpleInterval paddedLoc = region.getPaddedSpan();
            final Haplotype refHaplotype = AssemblyBasedCallerUtils.createReferenceHaplotype(region, paddedLoc, referenceReader);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, AssemblyBasedCallerUtils.createDummyStratifiedReadMap(refHaplotype, samplesList, readsHeader, region),
                    localGenotypingEngine.getPloidyModel(), Collections.emptyList(), hcArgs.standardArgs.genotypeArgs.supportVariants != null, VCpriors);
        }
        else {
            return NO_CALLS;
        }
    }

    /**
     * Shutdown this HC engine, closing resources as appropriate
     */
    public void shutdown() {
        minosSyncRandom();
        if (minosWarmActive) {
            minosTrailingRandom = minosLastRandom;
            minosTrailingWell = minosLastWell;
        }
        if (minosRepairOut != null) {
            minosRepairOut.close();
        }
        minosLikelihoodEngines.values().forEach(ReadLikelihoodCalculationEngine::close);
        for (final MinosGenotypingVariant v : minosGenotypingVariants) {
            if (v.writer != null) {
                v.writer.close();
            }
            try {
                if (v.diverged) {
                    Files.write(Paths.get(v.output + ".minos.regions"), v.regionLog);
                }
                Files.write(Paths.get(v.output + ".minos"), (v.diverged ? "diverged\n" : "exact\n").getBytes());
            } catch (final IOException e) {
                throw new UserException("cannot write status for " + v.output, e);
            }
        }
        minosWarmServe();
        likelihoodCalculationEngine.close();
        if (pdhmmLikelihoodCalculationEngine != null) pdhmmLikelihoodCalculationEngine.close();
        aligner.close();
        haplotypeBAMWriter.ifPresent(HaplotypeBAMWriter::close);
        assembledEventMapVcfOutputWriter.ifPresent(writer -> {assembledEventMapVariants.get().forEach(event -> writer.add(event)); writer.close();});
        if ( referenceReader != null) {
            try {
                referenceReader.close();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.close();
            } catch (IOException e) {
                throw new UserException("Error closing debug output stream", e);
            }
        }
        HaplotypeCallerGenotypingDebugger.close();
        // Write assembly region debug output if present
        assemblyEngine.printDebugHistograms();

    }

    protected Set<GATKRead> filterNonPassingReads( final AssemblyRegion activeRegion ) {
        // TODO: can we unify this additional filtering with makeStandardHCReadFilter()?

        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( AlignmentUtils.unclippedReadLength(rec) < READ_LENGTH_FILTER_THRESHOLD ||
                    rec.getMappingQuality() < hcArgs.mappingQualityThreshold ||
                    !ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(rec) ||
                    (hcArgs.keepRG != null && !rec.getReadGroup().equals(hcArgs.keepRG)) ) {
                if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                    HaplotypeCallerGenotypingDebugger.println("Filtered before assembly the read: " + rec);
                }
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    protected void removeReadsFromAllSamplesExcept(final String targetSample, final AssemblyRegion activeRegion) {
        final Set<GATKRead> readsToRemove = new LinkedHashSet<>();
        for( final GATKRead rec : activeRegion.getReads() ) {
            if( ! ReadUtils.getSampleName(rec, readsHeader).equals(targetSample) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos: batch evaluation of config variants in one run.
    //
    // -Dminos.gt.variants=<tsv>, one line per variant (evaluator/gatk.py variant_line): output.vcf.gz, then
    //   columns 1-4  genotyping: stand_call_conf, indel_heterozygosity, max_alternate_alleles, contamination
    //   columns 5-8  PairHMM: base_quality_score_threshold, pcr_indel_model, gap continuation, mismapping rate
    //   columns 9-14 assembly: min_pruning, min_dangling_branch_length, recover_all_dangling_branches,
    //                max_num_haplotypes_in_population, adaptive_pruning_initial_error_rate, pruning_lod_threshold
    //   columns 15-17 activity: active_probability_threshold, min/max_assembly_region_size
    // Each variant shares everything before the first stage it changes (genotyping engine; likelihoods; assembly;
    // region finding - see AssemblyRegionIterator). Each output must be bit-identical to a separate stock run of
    // that config; a variant whose random draws would differ is flagged diverged (and repaired by the runner).
    // ---------------------------------------------------------------------------------------------------------------
    private static final class MinosGenotypingVariant {
        final String output;
        final HaplotypeCallerGenotypingEngine engine;
        final String likelihoodKey;          // null = same PairHMM parameters as the main config
        final String assemblyKey;            // null = same assembly parameters as the main config
        final HaplotypeCallerArgumentCollection args;
        double[] activity = null;            // {threshold, min size, max size}; null = the main config's
        int stream = -1;                     // activity stream, -1 = the main config's regions
        final org.broadinstitute.hellbender.utils.MinosRecordingRandom rng;
        final org.broadinstitute.hellbender.utils.MinosRecordingWell well;
        final org.apache.commons.math3.random.RandomDataGenerator rdg;
        VariantContextWriter writer;
        boolean diverged = false;
        final List<String> regionLog = new ArrayList<>();

        /** One line per region: key, then this variant's draw positions before it and its draws in it. */
        void logRegion(final String key, final long[] before) {
            regionLog.add(key + "\t" + before[0] + "\t" + (rng.gaussians() - before[0]) + "\t" + before[1] + "\t"
                    + (well.steps() - before[1]));
        }

        MinosGenotypingVariant(final String output, final HaplotypeCallerGenotypingEngine engine, final String likelihoodKey,
                               final String assemblyKey, final HaplotypeCallerArgumentCollection args) {
            this(output, engine, likelihoodKey, assemblyKey, args,
                    Utils.minosGlobalRandom().copyState(), Utils.minosGlobalWell().copyState());
        }

        /** rng/well are the generator state at the start of the run, the same copy a variant built during init gets. */
        MinosGenotypingVariant(final String output, final HaplotypeCallerGenotypingEngine engine, final String likelihoodKey,
                               final String assemblyKey, final HaplotypeCallerArgumentCollection args,
                               final org.broadinstitute.hellbender.utils.MinosRecordingRandom rng,
                               final org.broadinstitute.hellbender.utils.MinosRecordingWell well) {
            this.args = args;
            this.output = output;
            this.engine = engine;
            this.likelihoodKey = likelihoodKey;
            this.assemblyKey = assemblyKey;
            this.rng = rng;
            this.well = well;
            this.rdg = new org.apache.commons.math3.random.RandomDataGenerator(well);
        }
    }

    private final List<MinosGenotypingVariant> minosGenotypingVariants = new ArrayList<>();
    private final Map<String, ReadLikelihoodCalculationEngine> minosLikelihoodEngines = new LinkedHashMap<>();
    private final Map<String, ReadThreadingAssembler> minosAssemblyEngines = new LinkedHashMap<>();

    private void minosInitGenotypingVariants() {
        final String spec = System.getProperty("minos.gt.variants");
        if (spec == null) {
            return;
        }
        if (!allCustomPloidies.isEmpty() || emitReferenceConfidence() || hcArgs.pileupDetectionArgs.generatePDHaplotypes) {
            throw new UserException("minos.gt.variants supports plain single-ploidy VCF mode only");
        }
        if (hcArgs.stepwiseFiltering || hcArgs.filterAlleles || hcArgs.dontGenotype
                || hcArgs.pileupDetectionArgs.usePileupDetection
                || assemblyDebugOutStream != null || hcArgs.justDetermineActiveRegions || hcArgs.sampleNameToUse != null
                || hcArgs.likelihoodArgs.likelihoodEngineImplementation != ReadLikelihoodCalculationEngine.Implementation.PairHMM) {
            throw new UserException("minos.gt.variants: a HaplotypeCaller mode the batch path does not replicate is enabled");
        }
        final List<String> lines;
        try {
            lines = Files.readAllLines(Paths.get(spec));
        } catch (final IOException e) {
            throw new UserException("cannot read " + spec, e);
        }
        for (final String line : lines) {
            if (line.isBlank() || line.startsWith("#")) {
                continue;
            }
            final String[] f = line.split("\t");
            final HaplotypeCallerArgumentCollection args = minosCloneArgs(hcArgs);
            args.standardArgs.genotypeArgs.standardConfidenceForCalling = Double.parseDouble(f[1]);
            args.standardArgs.genotypeArgs.indelHeterozygosity = Double.parseDouble(f[2]);
            args.standardArgs.genotypeArgs.maxAlternateAlleles = Integer.parseInt(f[3]);
            args.standardArgs.CONTAMINATION_FRACTION = Double.parseDouble(f[4]);
            args.standardArgs.setSampleContamination(java.util.Collections.emptyMap());   // rebuild as a fresh run would
            String likelihoodKey = null;
            if (f.length > 5) {                    // optional PairHMM columns: bqt, pcr_indel_model, gcp, mismapping
                args.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD = Byte.parseByte(f[5]);
                args.likelihoodArgs.pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCRErrorModel.valueOf(f[6]);
                args.likelihoodArgs.gcpHMM = Integer.parseInt(f[7]);
                args.likelihoodArgs.phredScaledGlobalReadMismappingRate = Integer.parseInt(f[8]);
                final String key = f[5] + "|" + f[6] + "|" + f[7] + "|" + f[8];
                final String mainKey = hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD + "|" + hcArgs.likelihoodArgs.pcrErrorModel
                        + "|" + hcArgs.likelihoodArgs.gcpHMM + "|" + hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate;
                if (!key.equals(mainKey)) {
                    likelihoodKey = key;
                    minosLikelihoodEngines.computeIfAbsent(key, k -> AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(
                            args.likelihoodArgs, args.fbargs, !args.softClipLowQualityEnds, args.likelihoodArgs.likelihoodEngineImplementation));
                }
            }
            String assemblyKey = null;
            if (f.length > 9) {                    // optional assembly columns
                args.assemblerArgs.minPruneFactor = Integer.parseInt(f[9]);
                args.assemblerArgs.minDanglingBranchLength = Integer.parseInt(f[10]);
                args.assemblerArgs.recoverAllDanglingBranches = Boolean.parseBoolean(f[11]);
                args.assemblerArgs.maxNumHaplotypesInPopulation = Integer.parseInt(f[12]);
                args.assemblerArgs.initialErrorRateForPruning = Double.parseDouble(f[13]);
                args.assemblerArgs.pruningLogOddsThreshold = Double.parseDouble(f[14]);
                final String key = minosAssemblyKey(args);
                if (!key.equals(minosAssemblyKey(hcArgs))) {
                    assemblyKey = key;
                    minosAssemblyEngines.computeIfAbsent(key, k -> args.createReadThreadingAssembler());
                }
            }
            final HaplotypeCallerGenotypingEngine engine =
                    new HaplotypeCallerGenotypingEngine(args, samplesList, !hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
            engine.setAnnotationEngine(annotationEngine);
            final MinosGenotypingVariant variant = new MinosGenotypingVariant(f[0], engine, likelihoodKey, assemblyKey, args);
            if (f.length > 15) {                   // optional activity columns: threshold, min size, max size
                variant.activity = new double[]{Double.parseDouble(f[15]), Integer.parseInt(f[16]), Integer.parseInt(f[17])};
                AssemblyRegionIterator.minosHook = this;
            }
            minosGenotypingVariants.add(variant);
        }
        logger.info("minos: " + minosGenotypingVariants.size() + " variants, " + minosLikelihoodEngines.size()
                + " extra likelihood engines");
        if (!minosGenotypingVariants.isEmpty()) {
            Utils.minosGlobalRandom().setRecording(true);
            Utils.minosGlobalWell().setRecording(true);
        }
    }

    /**
     * Random draws made outside genotyping (read downsampling during traversal) are shared by all variants in this
     * run. In each variant's own run they would come from that variant's stream, whose position differs whenever
     * its genotyping drew a different number of QD-jitter values. Replay them on each variant's generator: if any
     * value differs, that variant's own run would have downsampled differently, so its output here is not exact.
     *
     * A warm walk has no batch variants, but a later request replays these same draws (minosWarmOne). Still take
     * and keep them: otherwise every region stores an empty list, the request skips the downsampling draws, and
     * it reports exact while QD and QUAL differ from a fresh run.
     */
    private void minosSyncRandom() {
        if (minosGenotypingVariants.isEmpty() && !minosWarmActive) {
            return;
        }
        final List<int[]> draws = Utils.minosGlobalRandom().takeDraws();
        final List<int[]> wellDraws = Utils.minosGlobalWell().takeDraws();
        minosLastRandom = draws;
        minosLastWell = wellDraws;
        for (final MinosGenotypingVariant v : minosGenotypingVariants) {
            for (final int[] d : draws) {
                if (!v.rng.replay(d)) {
                    v.diverged = true;
                }
            }
            for (final int[] d : wellDraws) {
                if (!v.well.replay(d)) {
                    v.diverged = true;
                }
            }
        }
    }

    private void minosBeforeGenotyping() {
        if (minosGenotypingVariants.isEmpty() && !minosWarmActive) {
            return;
        }
        minosSyncRandom();
        Utils.minosGlobalRandom().setRecording(false);   // the main config's own genotyping draws are not shared
        Utils.minosGlobalWell().setRecording(false);
    }

    private void minosAfterGenotyping() {
        if (minosGenotypingVariants.isEmpty() && !minosWarmActive) {
            return;
        }
        Utils.minosGlobalRandom().setRecording(true);
        Utils.minosGlobalWell().setRecording(true);
    }

    /**
     * Deep copy of the argument collection, then a structural self-check against the original.
     * Java serialization cannot be used: the parent class AssemblyBasedCallerArgumentCollection is not Serializable
     * (its fields - minBaseQualityScore, assemblerArgs, likelihoodArgs, ... - would silently reset to defaults) and
     * LikelihoodEngineArgumentCollection holds a non-serializable PairHMMNativeArgumentCollection. GATK objects are
     * rebuilt field by field through their no-arg constructors; containers get new copies; everything else
     * (immutable values, feature inputs, paths) is shared, and never modified by the variant code.
     */
    static HaplotypeCallerArgumentCollection minosCloneArgs(final HaplotypeCallerArgumentCollection original) {
        final HaplotypeCallerArgumentCollection copy =
                (HaplotypeCallerArgumentCollection) minosDeepCopy(original, new java.util.IdentityHashMap<>());
        if (!minosDeepEquals(original, copy, new java.util.IdentityHashMap<>())) {
            throw new UserException("minos: HaplotypeCaller argument copy is not identical to the original");
        }
        return copy;
    }

    private static boolean minosIsValue(final Object o) {
        return o == null || o instanceof String || o instanceof Number || o instanceof Boolean || o instanceof Character
                || o instanceof Enum || o instanceof Class;
    }

    private static boolean minosIsGatkObject(final Object o) {
        return o.getClass().getName().startsWith("org.broadinstitute.hellbender.")
                && !(o instanceof org.broadinstitute.hellbender.engine.FeatureInput);
    }

    private static java.util.List<java.lang.reflect.Field> minosFields(final Class<?> cls) {
        final java.util.List<java.lang.reflect.Field> out = new ArrayList<>();
        for (Class<?> c = cls; c != null && c != Object.class; c = c.getSuperclass()) {
            for (final java.lang.reflect.Field f : c.getDeclaredFields()) {
                if (!java.lang.reflect.Modifier.isStatic(f.getModifiers())) {
                    f.setAccessible(true);
                    out.add(f);
                }
            }
        }
        return out;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static Object minosDeepCopy(final Object o, final java.util.IdentityHashMap<Object, Object> seen) {
        if (minosIsValue(o)) {
            return o;
        }
        if (seen.containsKey(o)) {
            return seen.get(o);
        }
        try {
            final Class<?> cls = o.getClass();
            if (cls.isArray()) {
                final int n = java.lang.reflect.Array.getLength(o);
                final Object a = java.lang.reflect.Array.newInstance(cls.getComponentType(), n);
                seen.put(o, a);
                for (int i = 0; i < n; i++) {
                    java.lang.reflect.Array.set(a, i, minosDeepCopy(java.lang.reflect.Array.get(o, i), seen));
                }
                return a;
            }
            if (o instanceof java.util.List) {
                final java.util.List l = new ArrayList();
                seen.put(o, l);
                for (final Object e : (java.util.List) o) {
                    l.add(minosDeepCopy(e, seen));
                }
                return l;
            }
            if (o instanceof java.util.Set) {
                final java.util.Set st = o instanceof java.util.SortedSet
                        ? new TreeSet(((java.util.SortedSet) o).comparator()) : new LinkedHashSet();
                seen.put(o, st);
                for (final Object e : (java.util.Set) o) {
                    st.add(minosDeepCopy(e, seen));
                }
                return st;
            }
            if (o instanceof java.util.Map) {
                final java.util.Map m = new LinkedHashMap();
                seen.put(o, m);
                for (final Object e : ((java.util.Map) o).entrySet()) {
                    final java.util.Map.Entry en = (java.util.Map.Entry) e;
                    m.put(minosDeepCopy(en.getKey(), seen), minosDeepCopy(en.getValue(), seen));
                }
                return m;
            }
            if (!minosIsGatkObject(o)) {
                return o;                              // shared (paths, feature inputs, ...): never modified
            }
            final java.lang.reflect.Constructor<?> ctor = cls.getDeclaredConstructor();
            ctor.setAccessible(true);
            final Object c = ctor.newInstance();
            seen.put(o, c);
            for (final java.lang.reflect.Field f : minosFields(cls)) {
                f.set(c, minosDeepCopy(f.get(o), seen));
            }
            return c;
        } catch (final ReflectiveOperationException e) {
            throw new UserException("minos: cannot copy " + o.getClass().getName(), e);
        }
    }

    @SuppressWarnings("rawtypes")
    private static boolean minosDeepEquals(final Object a, final Object b, final java.util.IdentityHashMap<Object, Object> seen) {
        if (minosIsValue(a) || minosIsValue(b)) {
            return java.util.Objects.equals(a, b);
        }
        if (seen.containsKey(a)) {
            return seen.get(a) == b;
        }
        seen.put(a, b);
        if (a.getClass() != b.getClass() && !(a instanceof java.util.Collection && b instanceof java.util.Collection)
                && !(a instanceof java.util.Map && b instanceof java.util.Map)) {
            return false;
        }
        try {
            if (a.getClass().isArray()) {
                final int n = java.lang.reflect.Array.getLength(a);
                if (n != java.lang.reflect.Array.getLength(b)) {
                    return false;
                }
                for (int i = 0; i < n; i++) {
                    if (!minosDeepEquals(java.lang.reflect.Array.get(a, i), java.lang.reflect.Array.get(b, i), seen)) {
                        return false;
                    }
                }
                return true;
            }
            if (a instanceof java.util.Collection) {
                final java.util.Iterator ia = ((java.util.Collection) a).iterator(), ib = ((java.util.Collection) b).iterator();
                while (ia.hasNext() && ib.hasNext()) {
                    if (!minosDeepEquals(ia.next(), ib.next(), seen)) {
                        return false;
                    }
                }
                return !ia.hasNext() && !ib.hasNext();
            }
            if (a instanceof java.util.Map) {
                return minosDeepEquals(new ArrayList(((java.util.Map) a).entrySet()).stream().map(e -> ((java.util.Map.Entry) e).getKey()).toArray(),
                        new ArrayList(((java.util.Map) b).entrySet()).stream().map(e -> ((java.util.Map.Entry) e).getKey()).toArray(), seen)
                    && minosDeepEquals(((java.util.Map) a).values().toArray(), ((java.util.Map) b).values().toArray(), seen);
            }
            if (!minosIsGatkObject(a)) {
                return a == b;
            }
            for (final java.lang.reflect.Field f : minosFields(a.getClass())) {
                if (!minosDeepEquals(f.get(a), f.get(b), seen)) {
                    return false;
                }
            }
            return true;
        } catch (final IllegalAccessException e) {
            throw new UserException("minos: cannot compare arguments", e);
        }
    }

    private static String minosAssemblyKey(final HaplotypeCallerArgumentCollection a) {
        return a.assemblerArgs.minPruneFactor + "|" + a.assemblerArgs.minDanglingBranchLength + "|"
                + a.assemblerArgs.recoverAllDanglingBranches + "|" + a.assemblerArgs.maxNumHaplotypesInPopulation + "|"
                + a.assemblerArgs.initialErrorRateForPruning + "|" + a.assemblerArgs.pruningLogOddsThreshold;
    }

    /**
     * Variants with their own assembly parameters. Called right after the main assembly, when `region` is
     * finalized (finalization depends only on traversal-level settings, and assembleReads skips it for a finalized
     * region) and not yet touched by anything else (trimming works on copies). For each distinct assembly setting:
     * assemble, then replicate the main path step by step (pileup events, trimming, read filtering, likelihoods per
     * PairHMM setting, realignment to the best haplotype, genotyping per variant). All random draws of a variant go
     * to its own generators; traversal draws so far are synced into them first.
     */
    private void minosCallAssemblyVariants(final AssemblyRegion region, final FeatureContext features,
                                           final ReferenceContext referenceContext, final List<Event> givenAlleles) {
        if (minosAssemblyEngines.isEmpty()) {
            return;
        }
        if (haplotypeBAMWriter.isPresent()) {      // created after initialize(), so checked here
            throw new UserException("minos.gt.variants: --bam-output is not supported with assembly variants");
        }
        minosSyncRandom();
        for (final Map.Entry<String, ReadThreadingAssembler> group : minosAssemblyEngines.entrySet()) {
            final List<MinosGenotypingVariant> members = new ArrayList<>();
            for (final MinosGenotypingVariant v : minosGenotypingVariants) {
                if (group.getKey().equals(v.assemblyKey) && v.stream < 0) {
                    members.add(v);
                }
            }
            if (members.isEmpty()) {
                continue;                          // all its variants have their own activity streams
            }
            // assembly and later steps draw no random numbers, but route any to the group's first variant anyway
            Utils.minosSetRandomOverride(members.get(0).rng);
            Utils.minosSetRandomDataOverride(members.get(0).rdg);
            final boolean recording = true;
            Utils.minosGlobalRandom().setRecording(false);
            Utils.minosGlobalWell().setRecording(false);
            try {
                minosCallAssemblyGroup(region, features, referenceContext, givenAlleles, group.getValue(), members);
            } finally {
                Utils.minosSetRandomOverride(null);
                Utils.minosSetRandomDataOverride(null);
                Utils.minosGlobalRandom().setRecording(recording);
                Utils.minosGlobalWell().setRecording(recording);
            }
        }
    }

    private void minosCallAssemblyGroup(final AssemblyRegion region, final FeatureContext features,
                                        final ReferenceContext referenceContext, final List<Event> givenAlleles,
                                        final ReadThreadingAssembler groupAssembler,
                                        final List<MinosGenotypingVariant> members) {
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(region, hcArgs, readsHeader,
                samplesList, logger, referenceReader, groupAssembler, aligner, !hcArgs.doNotCorrectOverlappingBaseQualities,
                hcArgs.fbargs, false);
        final Pair<Set<Event>, Set<Event>> goodAndBadPileupEvents =
                PileupBasedAlleles.goodAndBadPileupEvents(region.getAlignmentData(), hcArgs.pileupDetectionArgs, readsHeader, hcArgs.minBaseQualityScore);
        final Set<Event> goodPileupEvents = goodAndBadPileupEvents.getLeft();
        final Set<Event> badPileupEvents = goodAndBadPileupEvents.getRight();
        final SortedSet<Event> allVariationEvents = untrimmedAssemblyResult.getVariationEvents(hcArgs.maxMnpDistance).stream()
                .filter(event -> !badPileupEvents.contains(event))
                .collect(Collectors.toCollection(() -> new TreeSet<>(AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR)));
        goodPileupEvents.forEach(allVariationEvents::add);
        givenAlleles.forEach(allVariationEvents::add);
        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(region, allVariationEvents, referenceContext);
        if (!trimmingResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return;
        }
        final AssemblyResultSet assemblyResult = untrimmedAssemblyResult.trimTo(trimmingResult.getVariantRegion());
        assemblyResult.addGivenAlleles(givenAlleles, hcArgs.maxMnpDistance, aligner, hcArgs.getHaplotypeToReferenceSWParameters());
        assemblyResult.removeHaplotypesWithBadAlleles(hcArgs, badPileupEvents);
        assemblyResult.injectPileupEvents(region, hcArgs, aligner, goodPileupEvents);
        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        final List<GATKRead> readStubs = regionForGenotyping.getReads().stream()
                .filter(r -> AlignmentUtils.unclippedReadLength(r) < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        regionForGenotyping.removeAll(readStubs);
        final Collection<GATKRead> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKRead>> perSampleFilteredReadList = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, filteredReads);
        if (!assemblyResult.isVariationPresent() && !hcArgs.disableOptimizations) {
            return;
        }
        if (regionForGenotyping.size() == 0 && !hcArgs.disableOptimizations) {
            return;
        }
        final Map<String, List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, regionForGenotyping.getReads());
        final SWParameters readToHaplotypeSWParameters = hcArgs.getReadToHaplotypeSWParameters();
        final Map<String, AlleleLikelihoods<GATKRead, Haplotype>> byLikelihoodKey = new HashMap<>();
        for (final MinosGenotypingVariant v : members) {
            final String lk = v.likelihoodKey == null ? "" : v.likelihoodKey;
            AlleleLikelihoods<GATKRead, Haplotype> rl = byLikelihoodKey.get(lk);
            if (rl == null) {
                final ReadLikelihoodCalculationEngine engine = v.likelihoodKey == null
                        ? likelihoodCalculationEngine : minosLikelihoodEngines.get(v.likelihoodKey);
                rl = possiblyUncollapseHaplotypesInReadLikelihoods(untrimmedAssemblyResult,
                        engine.computeReadLikelihoods(assemblyResult, samplesList, reads, true));
                rl.changeEvidence(AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(rl,
                        assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner,
                        readToHaplotypeSWParameters));
                byLikelihoodKey.put(lk, rl);
            }
            Utils.minosSetRandomOverride(v.rng);
            Utils.minosSetRandomDataOverride(v.rdg);
            v.engine.assignGenotypeLikelihoods(rl.alleles(), rl, perSampleFilteredReadList,
                    assemblyResult.getFullReferenceWithPadding(), assemblyResult.getPaddedReferenceLoc(),
                    regionForGenotyping.getSpan(), features, givenAlleles, emitReferenceConfidence(), hcArgs.maxMnpDistance,
                    readsHeader, false, new HashSet<>(), rl)
                .getCalls().stream().map(RMSMappingQuality.getInstance()::finalizeRawMQ).forEach(v.writer::add);
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos: activity streams (see AssemblyRegionIterator). A stream's regions come straight from the traversal, so
    // each is called from scratch like callRegion does: the same early exits, then assembly and everything after it
    // per distinct assembly setting of the stream's variants.
    // ---------------------------------------------------------------------------------------------------------------
    private double[] minosMainActivity = null;
    private final List<List<MinosGenotypingVariant>> minosStreamMembers = new ArrayList<>();
    private final Map<String, ReadThreadingAssembler> minosStreamAssemblers = new LinkedHashMap<>();

    @Override
    public List<double[]> minosActivityStreams(final double activeProbThreshold, final int minRegionSize, final int maxRegionSize) {
        final double[] main = {activeProbThreshold, minRegionSize, maxRegionSize};
        if (minosMainActivity != null) {
            if (!Arrays.equals(minosMainActivity, main)) {
                throw new UserException("minos: activity settings changed between shards");
            }
        } else {
            minosMainActivity = main;
            final Map<String, Integer> index = new LinkedHashMap<>();
            for (final MinosGenotypingVariant v : minosGenotypingVariants) {
                if (v.activity == null || Arrays.equals(v.activity, main)) {
                    continue;
                }
                final String key = v.activity[0] + "|" + v.activity[1] + "|" + v.activity[2];
                v.stream = index.computeIfAbsent(key, k -> {
                    minosStreamMembers.add(new ArrayList<>());
                    return minosStreamMembers.size() - 1;
                });
                minosStreamMembers.get(v.stream).add(v);
                minosStreamAssemblers.computeIfAbsent(minosAssemblyKey(v.args), k -> v.args.createReadThreadingAssembler());
            }
            logger.info("minos: " + minosStreamMembers.size() + " activity streams");
        }
        final List<double[]> specs = new ArrayList<>();
        for (final List<MinosGenotypingVariant> members : minosStreamMembers) {
            specs.add(members.get(0).activity);
        }
        return specs;
    }

    @Override
    public void minosProcessRegion(final int stream, final AssemblyRegion region, final ReferenceContext referenceContext,
                                   final FeatureContext features) {
        final String key = minosRegionKey(region);
        final List<MinosGenotypingVariant> members = minosStreamMembers.get(stream);
        final List<long[]> before = new ArrayList<>();
        for (final MinosGenotypingVariant v : members) {
            before.add(new long[]{v.rng.gaussians(), v.well.steps()});
        }
        try {
            minosProcessStreamRegion(stream, region, referenceContext, features);
        } finally {
            for (int i = 0; i < members.size(); i++) {
                members.get(i).logRegion(key, before.get(i));
            }
        }
    }

    private void minosProcessStreamRegion(final int stream, final AssemblyRegion region, final ReferenceContext referenceContext,
                                          final FeatureContext features) {
        if (haplotypeBAMWriter.isPresent()) {
            throw new UserException("minos.gt.variants: --bam-output is not supported with activity variants");
        }
        if (!region.isActive()) {
            return;
        }
        final List<Event> givenAlleles = features.getValues(hcArgs.alleles).stream()
                .filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())
                .flatMap(vc -> GATKVariantContextUtils.splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
                .collect(Collectors.toList());
        if (givenAlleles.isEmpty() && region.size() == 0) {
            return;
        }
        minosSyncRandom();
        final Map<String, List<MinosGenotypingVariant>> groups = new LinkedHashMap<>();
        for (final MinosGenotypingVariant v : minosStreamMembers.get(stream)) {
            groups.computeIfAbsent(minosAssemblyKey(v.args), k -> new ArrayList<>()).add(v);
        }
        for (final Map.Entry<String, List<MinosGenotypingVariant>> group : groups.entrySet()) {
            final List<MinosGenotypingVariant> members = group.getValue();
            Utils.minosSetRandomOverride(members.get(0).rng);
            Utils.minosSetRandomDataOverride(members.get(0).rdg);
            Utils.minosGlobalRandom().setRecording(false);
            Utils.minosGlobalWell().setRecording(false);
            try {
                minosCallAssemblyGroup(region, features, referenceContext, givenAlleles,
                        minosStreamAssemblers.get(group.getKey()), members);
            } finally {
                Utils.minosSetRandomOverride(null);
                Utils.minosSetRandomDataOverride(null);
                Utils.minosGlobalRandom().setRecording(true);
                Utils.minosGlobalWell().setRecording(true);
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos: repair of a diverged batch variant. The variant's own run (own random draws, so its read pass keeps exactly
    // the reads its real run keeps) with -Dminos.repair.log=<the batch's region log for it>. A region whose span,
    // activity and reads equal a logged region gives the same calls as in the batch, so it is skipped: no assembly,
    // likelihoods or genotyping, the generators are advanced by the draws the batch logged for it. That is only valid
    // while the generators sit exactly where the batch's did before the region; the first time they do not, every
    // later region is computed. -Dminos.repair.out=<file> lists each region as S(kipped) or C(omputed) for the splice.
    // ---------------------------------------------------------------------------------------------------------------
    private Map<String, long[]> minosRepairLog = null;
    private java.io.PrintWriter minosRepairOut = null;
    private boolean minosRepairComputeAll = false;

    private void minosInitRepair() {
        final String spec = System.getProperty("minos.repair.log");
        if (spec == null) {
            return;
        }
        if (System.getProperty("minos.gt.variants") != null) {
            throw new UserException("minos.repair.log and minos.gt.variants cannot be combined");
        }
        minosRepairLog = new HashMap<>();
        try {
            for (final String line : Files.readAllLines(Paths.get(spec))) {
                final String[] f = line.split("\t");
                minosRepairLog.put(f[0], new long[]{Long.parseLong(f[1]), Long.parseLong(f[2]), Long.parseLong(f[3]), Long.parseLong(f[4])});
            }
            minosRepairOut = new java.io.PrintWriter(Files.newBufferedWriter(Paths.get(System.getProperty("minos.repair.out"))));
        } catch (final IOException e) {
            throw new UserException("cannot read repair log " + spec, e);
        }
    }

    private List<VariantContext> minosRepairRegion(final AssemblyRegion region, final FeatureContext features,
                                                   final ReferenceContext referenceContext) {
        final String key = minosRegionKey(region);
        final long[] logged = minosRepairLog.get(key);
        if (logged != null && !minosRepairComputeAll) {
            if (Utils.minosGlobalRandom().gaussians() == logged[0] && Utils.minosGlobalWell().steps() == logged[2]) {
                Utils.minosGlobalRandom().skipGaussians(logged[1]);
                Utils.minosGlobalWell().skipSteps(logged[3]);
                minosRepairOut.println("S\t" + key);
                return NO_CALLS;
            }
            minosRepairComputeAll = true;        // draw positions differ from here on: compute every later region
        }
        minosRepairOut.println("C\t" + key);
        return minosCallRegion(region, features, referenceContext);
    }

    /** Span, activity and a fingerprint of the region's reads, taken before anything touches the region. */
    private static String minosRegionKey(final AssemblyRegion region) {
        long h = 1125899906842597L;
        for (final GATKRead r : region.getReads()) {
            h = 31 * h + r.getName().hashCode();
            h = 31 * h + r.getStart();
            h = 31 * h + r.getEnd();
            h = 31 * h + r.getCigar().toString().hashCode();
            h = 31 * h + r.getMappingQuality();
            h = 31 * h + (r.isReverseStrand() ? 1 : 0) + (r.isFirstOfPair() ? 2 : 0) + (r.isSecondOfPair() ? 4 : 0);
            h = 31 * h + Arrays.hashCode(r.getBasesNoCopy());
            h = 31 * h + Arrays.hashCode(r.getBaseQualitiesNoCopy());
        }
        return region.getSpan().getContig() + ":" + region.getSpan().getStart() + "-" + region.getSpan().getEnd() + "|"
                + region.isActive() + "|" + region.getReads().size() + "|" + Long.toHexString(h);
    }

    // ---------------------------------------------------------------------------------------------------------------
    // minos warm JVM: after the traversal, stay up and genotype later variants of this same run.
    // Only genotyping parameters may change. The stored likelihoods are the ones an in-batch variant would
    // have genotyped; a variant whose read-pass draws would differ is flagged diverged, same as a batch.
    // -Dminos.warm.dir=<dir written by the runner>. Ready file "ok" or "no". Requests are req-00001.tsv.
    // ---------------------------------------------------------------------------------------------------------------
    private java.nio.file.Path minosWarmDir;
    private boolean minosWarmActive;
    private String minosWarmKey;
    private org.broadinstitute.hellbender.utils.MinosRecordingRandom minosWarmBaseRandom;
    private org.broadinstitute.hellbender.utils.MinosRecordingWell minosWarmBaseWell;
    private List<int[]> minosLastRandom = java.util.Collections.emptyList();
    private List<int[]> minosLastWell = java.util.Collections.emptyList();
    private List<int[]> minosTrailingRandom = java.util.Collections.emptyList();
    private List<int[]> minosTrailingWell = java.util.Collections.emptyList();
    private final List<WarmRegion> minosWarmRegions = new ArrayList<>();
    private static final int MINOS_WARM_MAX_REGIONS = 200_000;   // one chr18 round is ~64k regions; the heap guard is the real limit
    private static final long MINOS_WARM_MIN_FREE = 256L * 1024 * 1024;
    private static final long MINOS_WARM_IDLE_MS = 15 * 60 * 1000L;

    private static final class WarmRegion {
        String key;
        List<int[]> randomDraws;
        List<int[]> wellDraws;
        AlleleLikelihoods<GATKRead, Haplotype> likelihoods;
        Map<String, List<GATKRead>> filteredReads;
        byte[] fullRef;
        SimpleInterval paddedLoc;
        SimpleInterval span;
        FeatureContext features;
        List<Event> givenAlleles;
        SAMFileHeader readsHeader;
        Set<Integer> suspicious;
    }

    private void minosInitWarm() {
        final String dir = System.getProperty("minos.warm.dir");
        if (dir == null) {
            return;
        }
        minosWarmDir = java.nio.file.Paths.get(dir);
        if (!allCustomPloidies.isEmpty() || emitReferenceConfidence() || hcArgs.pileupDetectionArgs.generatePDHaplotypes
                || hcArgs.stepwiseFiltering || hcArgs.filterAlleles || hcArgs.dontGenotype
                || hcArgs.pileupDetectionArgs.usePileupDetection || assemblyDebugOutStream != null
                || hcArgs.justDetermineActiveRegions || hcArgs.sampleNameToUse != null
                || hcArgs.likelihoodArgs.likelihoodEngineImplementation != ReadLikelihoodCalculationEngine.Implementation.PairHMM) {
            return;   // active stays false; shutdown writes "no" and exits
        }
        minosWarmActive = true;
        AssemblyRegionIterator.minosHook = this;   // records the main activity settings; adds no stream when none differ
        minosWarmBaseRandom = Utils.minosGlobalRandom().copyState();
        minosWarmBaseWell = Utils.minosGlobalWell().copyState();
        Utils.minosGlobalRandom().setRecording(true);
        Utils.minosGlobalWell().setRecording(true);
    }

    private void minosWarmCapture(final AlleleLikelihoods<GATKRead, Haplotype> likelihoods,
                                  final Map<String, List<GATKRead>> filteredReads, final byte[] fullRef,
                                  final SimpleInterval paddedLoc, final SimpleInterval span, final FeatureContext features,
                                  final List<Event> givenAlleles, final SAMFileHeader readsHeader,
                                  final Set<Integer> suspicious) {
        if (!minosWarmActive) {
            return;
        }
        final long free = Runtime.getRuntime().maxMemory() - Runtime.getRuntime().totalMemory()
                + Runtime.getRuntime().freeMemory();
        if (minosWarmRegions.size() >= MINOS_WARM_MAX_REGIONS || free < MINOS_WARM_MIN_FREE) {
            logger.warn("minos warm: dropping snapshots at " + minosWarmRegions.size()
                    + " regions, " + (free / (1024 * 1024)) + " MiB free");
            minosWarmActive = false;
            minosWarmRegions.clear();
            return;
        }
        final WarmRegion region = new WarmRegion();
        region.key = minosWarmKey;
        region.randomDraws = minosLastRandom;
        region.wellDraws = minosLastWell;
        region.likelihoods = likelihoods;
        region.filteredReads = filteredReads;
        region.fullRef = fullRef == null ? null : fullRef.clone();
        region.paddedLoc = paddedLoc;
        region.span = span;
        region.features = features;
        region.givenAlleles = givenAlleles;
        region.readsHeader = readsHeader;
        region.suspicious = new HashSet<>(suspicious);
        minosWarmRegions.add(region);
    }

    private void minosWarmServe() {
        if (minosWarmDir == null) {
            return;
        }
        try {
            Files.createDirectories(minosWarmDir);
            if (!minosWarmActive || minosWarmRegions.isEmpty() || minosSequenceDictionary == null) {
                logger.info("minos warm: not keeping this JVM, regions=" + minosWarmRegions.size());
                Files.writeString(minosWarmDir.resolve("ready"), "no\n");
                return;
            }
            logger.info("minos warm: keeping " + minosWarmRegions.size() + " regions");
            Files.writeString(minosWarmDir.resolve("ready"), "ok\n");
            minosWarmLastRequest = System.currentTimeMillis();
            for (int n = 1; ; n++) {
                final java.nio.file.Path req = minosWarmDir.resolve(String.format("req-%05d.tsv", n));
                while (!Files.exists(req)) {
                    if (Files.exists(minosWarmDir.resolve("stop"))) {
                        return;
                    }
                    if (minosWarmIdle(n)) {
                        return;
                    }
                    Thread.sleep(20);
                }
                final String line = Files.readString(req).trim();
                if (line.equals("BYE")) {
                    Files.writeString(minosWarmDir.resolve(String.format("req-%05d.done", n)), "bye\n");
                    return;
                }
                minosWarmOne(n, line);
            }
        } catch (final Exception e) {
            logger.warn("minos warm session ended: " + e.getMessage());
            try {
                Files.writeString(minosWarmDir.resolve("ready"), "no\n");
            } catch (final IOException ignored) {
                // the runner treats a dead process as a cold run
            }
        }
    }

    private long minosWarmLastRequest = System.currentTimeMillis();

    private boolean minosWarmIdle(final int n) {
        return System.currentTimeMillis() - minosWarmLastRequest > MINOS_WARM_IDLE_MS;
    }

    private void minosWarmOne(final int n, final String line) throws IOException {
        minosWarmLastRequest = System.currentTimeMillis();
        final java.nio.file.Path done = minosWarmDir.resolve(String.format("req-%05d.done", n));
        final String[] f = line.split("\t");
        if (f.length < 18) {
            Files.writeString(done, "refused\n");
            return;
        }
        final String mainPair = hcArgs.likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD + "|" + hcArgs.likelihoodArgs.pcrErrorModel
                + "|" + hcArgs.likelihoodArgs.gcpHMM + "|" + hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate;
        final String pair = f[5] + "|" + f[6] + "|" + f[7] + "|" + f[8];
        final HaplotypeCallerArgumentCollection probe = minosCloneArgs(hcArgs);
        probe.assemblerArgs.minPruneFactor = Integer.parseInt(f[9]);
        probe.assemblerArgs.minDanglingBranchLength = Integer.parseInt(f[10]);
        probe.assemblerArgs.recoverAllDanglingBranches = Boolean.parseBoolean(f[11]);
        probe.assemblerArgs.maxNumHaplotypesInPopulation = Integer.parseInt(f[12]);
        probe.assemblerArgs.initialErrorRateForPruning = Double.parseDouble(f[13]);
        probe.assemblerArgs.pruningLogOddsThreshold = Double.parseDouble(f[14]);
        final double[] wantActivity = {Double.parseDouble(f[15]), Integer.parseInt(f[16]), Integer.parseInt(f[17])};
        final boolean sameActivity = minosMainActivity != null && Arrays.equals(minosMainActivity, wantActivity);
        if (!pair.equals(mainPair) || !minosAssemblyKey(probe).equals(minosAssemblyKey(hcArgs)) || !sameActivity) {
            Files.writeString(done, "refused\n");
            return;
        }
        final HaplotypeCallerArgumentCollection args = minosCloneArgs(hcArgs);
        args.standardArgs.genotypeArgs.standardConfidenceForCalling = Double.parseDouble(f[1]);
        args.standardArgs.genotypeArgs.indelHeterozygosity = Double.parseDouble(f[2]);
        args.standardArgs.genotypeArgs.maxAlternateAlleles = Integer.parseInt(f[3]);
        args.standardArgs.CONTAMINATION_FRACTION = Double.parseDouble(f[4]);
        args.standardArgs.setSampleContamination(java.util.Collections.emptyMap());
        final HaplotypeCallerGenotypingEngine engine =
                new HaplotypeCallerGenotypingEngine(args, samplesList, !hcArgs.doNotRunPhysicalPhasing, hcArgs.applyBQD);
        engine.setAnnotationEngine(annotationEngine);
        final MinosGenotypingVariant variant = new MinosGenotypingVariant(f[0], engine, null, null, args,
                minosWarmBaseRandom.copyState(), minosWarmBaseWell.copyState());
        variant.writer = GATKVariantContextUtils.createVCFWriter(java.nio.file.Paths.get(f[0]), minosSequenceDictionary,
                false, htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY);
        variant.writer.writeHeader(makeVCFHeader(minosSequenceDictionary, minosDefaultHeaderLines));
        try {
            for (final WarmRegion region : minosWarmRegions) {
                final long[] before = new long[]{variant.rng.gaussians(), variant.well.steps()};
                for (final int[] d : region.randomDraws) {
                    if (!variant.rng.replay(d)) {
                        variant.diverged = true;
                    }
                }
                for (final int[] d : region.wellDraws) {
                    if (!variant.well.replay(d)) {
                        variant.diverged = true;
                    }
                }
                Utils.minosSetRandomOverride(variant.rng);
                Utils.minosSetRandomDataOverride(variant.rdg);
                try {
                    engine.assignGenotypeLikelihoods(region.likelihoods.alleles(), region.likelihoods, region.filteredReads,
                            region.fullRef, region.paddedLoc, region.span, region.features, region.givenAlleles,
                            false, hcArgs.maxMnpDistance, region.readsHeader, false, region.suspicious, region.likelihoods)
                            .getCalls().stream().map(RMSMappingQuality.getInstance()::finalizeRawMQ).forEach(variant.writer::add);
                } finally {
                    Utils.minosSetRandomOverride(null);
                    Utils.minosSetRandomDataOverride(null);
                }
                variant.logRegion(region.key, before);
            }
            for (final int[] d : minosTrailingRandom) {
                if (!variant.rng.replay(d)) {
                    variant.diverged = true;
                }
            }
            for (final int[] d : minosTrailingWell) {
                if (!variant.well.replay(d)) {
                    variant.diverged = true;
                }
            }
        } finally {
            variant.writer.close();
        }
        if (variant.diverged) {
            Files.write(java.nio.file.Paths.get(variant.output + ".minos.regions"), variant.regionLog);
        }
        Files.writeString(java.nio.file.Paths.get(variant.output + ".minos"), (variant.diverged ? "diverged\n" : "exact\n"));
        Files.writeString(done, (variant.diverged ? "diverged\n" : "exact\n"));
    }
}
