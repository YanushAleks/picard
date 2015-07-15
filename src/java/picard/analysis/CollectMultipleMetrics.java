package picard.analysis;

import htsjdk.samtools.util.CollectionUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * Class that is designed to instantiate and execute multiple metrics programs that extend
 * SinglePassSamProgram while making only a single pass through the SAM file and supplying
 * each program with the records as it goes.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Takes an input BAM and reference sequence and runs one or more Picard " +
                "metrics modules at the same time to cut down on I/O.  These include CollectAlignmentSummaryMetrics, " +
                "CollectInsertSizeMetrics, QualityScoreDistribution,  MeanQualityByCycle, and " +
                "CollectBaseDistributionByCycle.  Produces a pdf and a txt file for each tool with the exception of the " +
                "CollectAlignmentSummaryMetrics tool, which outputs only a txt file.  To view the \"txt\" files, add " +
                " \".txt\" to the suffix of each output file (but not the \".pdf\" files)." +
                "<br /><br />" +
                "" +
                "Currently all programs are run with default options and fixed output extensions, " +
                "but this may become more flexible in future.  Reference sequence file is required." +
                "<br />" +
                "<h4>Usage example:</h4>" +
                "<pre>" +
                "java -jar picard.jar MeanQualityByCycle \\<br />" +
                "     -I=input.bam \\<br />" +
                "     -O=MM_output.txt \\<br />" +
                "     -R=Reference.fasta <br /> " +
                "</pre>" +
                "<hr />"
        ,

        usageShort = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM",
        programGroup = Metrics.class
)
public class CollectMultipleMetrics extends CommandLineProgram {

    /**
     * This interface allows developers to create Programs to run in addition to the ones defined in the Program enum.
     */
    public static interface ProgramInterface {
        SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference);
    }

    public static enum Program implements ProgramInterface {
        CollectAlignmentSummaryMetrics {
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final CollectAlignmentSummaryMetrics program = new CollectAlignmentSummaryMetrics();
                program.OUTPUT = new File(outbase + ".alignment_summary_metrics");

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        },
        CollectInsertSizeMetrics {
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final CollectInsertSizeMetrics program = new CollectInsertSizeMetrics();
                program.OUTPUT = new File(outbase + ".insert_size_metrics");
                program.Histogram_FILE = new File(outbase + ".insert_size_histogram.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        },
        QualityScoreDistribution {
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final QualityScoreDistribution program = new QualityScoreDistribution();
                program.OUTPUT = new File(outbase + ".quality_distribution_metrics");
                program.CHART_OUTPUT = new File(outbase + ".quality_distribution.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        },
        MeanQualityByCycle {
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final MeanQualityByCycle program = new MeanQualityByCycle();
                program.OUTPUT = new File(outbase + ".quality_by_cycle_metrics");
                program.CHART_OUTPUT = new File(outbase + ".quality_by_cycle.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        },
        CollectBaseDistributionByCycle {
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final CollectBaseDistributionByCycle program = new CollectBaseDistributionByCycle();
                program.OUTPUT = new File(outbase + ".base_distribution_by_cycle_metrics");
                program.CHART_OUTPUT = new File(outbase + ".base_distribution_by_cycle.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        },
        CollectGcBiasMetrics {
            public SinglePassSamProgram makeInstance(final String outbase, final File input, final File reference) {
                final CollectGcBiasMetrics program = new CollectGcBiasMetrics();
                program.OUTPUT = new File(outbase + ".gc_bias.detail_metrics");
                program.SUMMARY_OUTPUT = new File(outbase + ".gc_bias.summary_metrics");
                program.CHART_OUTPUT = new File(outbase + ".gc_bias.pdf");

                program.INPUT = input;
                program.METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS,
                        MetricAccumulationLevel.LIBRARY);
                program.WINDOW_SIZE = 100;
                program.MINIMUM_GENOME_FRACTION = 1.0E-5;
                program.IS_BISULFITE_SEQUENCED = false;
                program.ASSUME_SORTED = false;

                //GC_Bias actually uses the class-level REFERENCE_SEQUENCE variable.
                program.REFERENCE_SEQUENCE = reference;

                return program;
            }
        }

    }

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;


    @Option(doc = "If true (default), then the sort order in the header file will be ignored.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = true;

    @Option(doc = "Stop after processing N reads, mainly for debugging.")
    public int STOP_AFTER = 0;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Base name of output files.")
    public String OUTPUT;

    @Option(doc = "List of metrics programs to apply during the pass through the SAM file.")
    public List<Program> PROGRAM = CollectionUtil.makeList(Program.values());

    /**
     * Contents of PROGRAM list is transferred to this list during command-line validation, so that an outside
     * developer can invoke this class programmatically and provide alternative Programs to run by calling
     * setProgramsToRun().
     */
    private List<ProgramInterface> programsToRun;

    // Stock main method
    public static void main(final String[] args) {
        new CollectMultipleMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        programsToRun = new ArrayList<ProgramInterface>(PROGRAM);
        return super.customCommandLineValidation();
    }

    /**
     * Use this method when invoking CollectMultipleMetrics programmatically to run programs other than the ones
     * available via enum.  This must be called before doWork().
     */
    public void setProgramsToRun(List<ProgramInterface> programsToRun) {
        this.programsToRun = programsToRun;
    }

    @Override
    public int doWork() {
        if (OUTPUT.endsWith(".")) {
            OUTPUT = OUTPUT.substring(0, OUTPUT.length() - 1);
        }

        final List<SinglePassSamProgram> programs = new ArrayList<SinglePassSamProgram>();
        for (final ProgramInterface program : new HashSet<ProgramInterface>(programsToRun)) {
            final SinglePassSamProgram instance = program.makeInstance(OUTPUT , INPUT, REFERENCE_SEQUENCE);

            // Generally programs should not be accessing these directly but it might make things smoother
            // to just set them anyway
            instance.INPUT = INPUT;
            instance.REFERENCE_SEQUENCE = REFERENCE_SEQUENCE;

            instance.setDefaultHeaders(getDefaultHeaders());

            programs.add(instance);
        }

        SinglePassSamProgram.makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, programs);

        return 0;
    }
}