from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars
from os import cpu_count, system

from bio_types import *
from bio_methods import *
from common import pipelineDriver

fallback_warning_shown = False


def getFileNames(options: OptionsDict) -> FastqSet:
    """
    Returns the input / source filename from the options. First, attempts the --fastq
    option as a full path. If that doesn't exist, construct a path from --fastq-dir
    as well and check it. If neither are present, this function exits the process.
    """
    global fallback_warning_shown

    fastq = options["fastq"]
    fastq_dir = options["fastq_dir"]

    if exists(expandvars(fastq)) == True:
        return [fastq]

    filename = "{FASTQ_DIR}/{FASTQ}".format(FASTQ_DIR=fastq_dir, FASTQ=fastq)

    if exists(expandvars(filename)) == True:
        return [filename]

    print("Check your --fastq and --fastq-dir parameters")
    quit(1)


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return ["{PIPELINE}/{SAMPLE}.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample)]


def runIdentityPreprocessor(
    script: TextIOWrapper,
    r1: str,
    o1: str,
):
    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} ]]; then
    logthis "${{yellow}}Running identity preprocessor${{reset}}"

    ln -s {R1} {O1}

    logthis "${{yellow}}Identity preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            O1=o1,
        )
    )


def runCutadaptPreprocessor(
    script: TextIOWrapper,
    r1: str,
    o1: str,
    options: OptionsDict,
):
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]

    script.write(
        """
#
# run the cutadapt preprocessor
#
if [[ ! -f {O1} ]]; then
    logthis "${{yellow}}Running cutadapt preprocessor${{reset}}"

    cutadapt \\
        --cores {THREADS} \\
        --length 300 \\
        --cut 15 \\
        --minimum-length 30 \\
        --quality-cutoff 18,20 \\
        --trim-n \\
        --adapter ATCACCGACTGCCCATAGAGAGGCTGAGAC --times 3 \\
        --adapter AAAAAAAAAA$ --times 3 \\
        --adapter GGGGGGGGGG$ --times 3 \\
        --output {O1} \\
        {R1} >{STATS}/{SAMPLE}.cutadapt.log

    logthis "${{yellow}}cutadapt preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1, O1=o1, SAMPLE=sample, STATS=stats, THREADS=threads
        )
    )


def preprocessFASTQ(
    script: TextIOWrapper,
    r1: str,
    o1: str,
    options: OptionsDict,
):
    preprocessor = options["preprocessor"]

    if preprocessor == "none":
        runIdentityPreprocessor(script, r1, o1)
    elif preprocessor == "cutadapt":
        runCutadaptPreprocessor(script, r1, o1, options)
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)


def runBwaAligner(
    script: TextIOWrapper,
    r1: str,
    options: OptionsDict,
):
    aligner = options["aligner"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    {ALIGNER} mem \\
        -R "@RG\\tID:{SAMPLE}\\tPL:IONTORRENT\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        -t {THREADS} \\
        -Y \\
        -M \\
        -v 1 \\
        -S \\
        -P \\
        {REFERENCE}/{ASSEMBLY}.fna \\
        {R1} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            ALIGNER=aligner,
            R1=r1,
            ASSEMBLY=assembly,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            BIN=bin,
        )
    )


def preprocessAndAlign(script: TextIOWrapper, options: OptionsDict):
    filenames = getFileNames(options)
    trimmedFilenames = getTrimmedFileNames(options)

    preprocessFASTQ(
        script,
        filenames[0],
        trimmedFilenames[0],
        options,
    )
    runBwaAligner(script, trimmedFilenames[0], options)


def defineArguments(panel_choices: List[str], panel_choice_help: str) -> Namespace:
    parser = ArgumentParser()
    parser.set_defaults(doQC=False, cleanIntermediateFiles=True)

    parser.add_argument(
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="Short name of sample, all files will use this as an alias",
    )
    parser.add_argument(
        "--fastq",
        required=True,
        action="store",
        metavar="FASTQ",
        dest="fastq",
        help="Path name of FASTQ relative to --fastq-dir",
    )

    parser.add_argument(
        "--work-dir",
        required=True,
        action="store",
        metavar="WORKING_DIR",
        dest="working",
        help="Working directory, e.g. base for $WORKING/pipeline, $WORKING/stats",
    )

    parser.add_argument(
        "--run-qc",
        action="store_true",
        dest="runQc",
        default=False,
        help="Enable running any / all QC processes",
    )
    parser.add_argument(
        "--skip-multi-qc",
        action="store_false",
        dest="doMultiQc",
        default=True,
        help="Skip running multiqc QC process on input and output files",
    )
    parser.add_argument(
        "--skip-variant-qc",
        action="store_false",
        dest="doVariantQc",
        default=True,
        help="Skip running variant QC process on input and output files",
    )
    parser.add_argument(
        "--skip-alignment-qc",
        action="store_false",
        dest="doAlignmentQc",
        default=True,
        help="Skip alignment QC process on input and output files",
    )
    parser.add_argument(
        "--skip-picard-qc",
        action="store_false",
        dest="doPicardQc",
        default=True,
        help="Skip any picard QC process on input and output files",
    )
    parser.add_argument(
        "--skip-fastqc",
        action="store_true",
        dest="skipFastQc",
        default=False,
        help="Skip running FASTQC statistics",
    )

    parser.add_argument(
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["cutadapt", "none"],
        help="Optionally run a FASTQ preprocessor",
    )

    parser.add_argument(
        "--unmapped",
        action="store_true",
        dest="processUnmapped",
        default=False,
        help="Extract unmapped reads into secondary _R1 and _R2 FASTQ files",
    )

    parser.add_argument(
        "--align-only",
        action="store_true",
        dest="alignOnly",
        default=False,
        help="Only run alignment and sorting processes",
    )

    parser.add_argument(
        "--aligner",
        action="store",
        dest="aligner",
        default="bwa",
        choices=["bwa", "bwa-mem2"],
        help="Use 'bwa' or 'bwa-mem2' as the aligner.",
    )

    parser.add_argument(
        "--sorter",
        action="store",
        dest="sorter",
        default="biobambam",
        choices=["biobambam", "samtools"],
        help="Use 'biobambam' or 'samtools' as the sorter.",
    )

    parser.add_argument(
        "--skip-annotation",
        action="store_true",
        dest="skipAnnotation",
        default=False,
        help="Skip all VCF annotation processes",
    )

    parser.add_argument(
        "--no-color",
        action="store_true",
        dest="noColor",
        default=False,
        help="Turn off colorized log output",
    )

    parser.add_argument(
        "--reference-assembly",
        required=True,
        action="store",
        metavar="REFERENCE_ASSEMBLY",
        dest="referenceAssembly",
        help="Base name of the reference assembly",
    )

    parser.add_argument(
        "--reference-name",
        required=True,
        action="store",
        metavar="REFERENCE_NAME",
        dest="referenceName",
        choices=panel_choices,
        help=panel_choice_help,
    )

    parser.add_argument(
        "--reference-dir",
        action="store",
        metavar="REFERENCE_DIR",
        dest="reference",
        help="Location of the reference genome files",
    )
    parser.add_argument(
        "--pipeline-dir",
        action="store",
        metavar="PIPELINE_DIR",
        dest="pipeline",
        help="Output directory for processing",
    )
    parser.add_argument(
        "--fastq-dir",
        action="store",
        metavar="FASTQ_DIR",
        dest="fastq_dir",
        help="Location of R1 and R2 files",
    )
    parser.add_argument(
        "--temp-dir",
        action="store",
        metavar="TEMP_DIR",
        dest="temp",
        help="Temporary storage for alignment",
    )
    parser.add_argument(
        "--stats-dir",
        action="store",
        metavar="STATS_DIR",
        dest="stats",
        help="Destination for statistics and QC files",
    )
    parser.add_argument(
        "--bin-dir",
        action="store",
        metavar="BIN_DIR",
        dest="bin",
        help="Install location of all tooling",
    )
    parser.add_argument(
        "--lib-dir",
        action="store",
        metavar="LIB_DIR",
        dest="lib",
        help="Install location of all libraries",
    )

    parser.add_argument(
        "--adapter-fasta",
        action="store",
        metavar="ADAPTER_FASTA",
        dest="adapters",
        default="",
        help="Optional FASTA file containing list of adapters to pass to fastp",
    )

    parser.add_argument(
        "--script",
        action="store",
        metavar="SHELL_SCRIPT",
        dest="script",
        default="pipeline-runner",
        help="Filename of bash shell to create",
    )

    parser.add_argument(
        "--cores",
        action="store",
        dest="cores",
        metavar="CPU_COUNT",
        default=cpu_count(),
        help="Specify the number of available CPU",
    )

    parser.add_argument(
        "--read-limit",
        action="store",
        dest="read-limit",
        metavar="READ_LIMIT",
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    return parser.parse_args()


def verifyFileNames(options: OptionsDict):
    fastq = options["fastq"]
    fastq_dir = options["fastq_dir"]

    if exists(expandvars(fastq)) == True:
        return fastq

    filename = "{FASTQ_DIR}/{FASTQ}".format(FASTQ_DIR=fastq_dir, FASTQ=fastq)

    if exists(expandvars(filename)) == True:
        return filename

    print("Check your --fastq and --fastq-dir parameters")
    quit(1)


def runSingleEndPipeline(script: TextIOWrapper, options: OptionsDict, filenames):
    commonPipeline(script, options, filenames, preprocessAndAlign)


def main(panel_choices: List[str], panel_choice_help: str):
    pipelineDriver(
        panel_choices,
        panel_choice_help,
        defineArguments,
        verifyFileNames,
        getFileNames,
        runSingleEndPipeline,
    )
