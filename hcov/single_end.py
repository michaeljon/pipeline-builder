from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars
from os import cpu_count, system

from bio_types import *
from common import *
from bio_methods import *

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

    return [
        "{PIPELINE}/{SAMPLE}.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample)
    ]


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
        print(
            "Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(
                PREPROCESSOR=preprocessor
            )
        )
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


def runHisatAligner(
    script: TextIOWrapper,
    r1: str,
    options: OptionsDict,
):
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

    hisat2 \\
            -x {REFERENCE}/{ASSEMBLY} \\
            -1 {R1} \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:IONTORRENT" \\
            --rg "PU:unspecified" \\
            --rg "LB:{SAMPLE}" \\
            --rg "SM:{SAMPLE}" \\
            --no-spliced-alignment \\
            --no-unal \\
            --threads {THREADS} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            BIN=bin,
        )
    )


def alignFASTQ(
    script: TextIOWrapper,
    r1: str,
    options: OptionsDict,
):
    aligner = options["aligner"]

    if aligner == "bwa" or aligner == "bwa-mem2":
        runBwaAligner(script, r1, options)
    elif aligner == "hisat2":
        runHisatAligner(script, r1, options)
    else:
        print(
            "Unexpected value {ALIGNER} given for the --aligner option".format(
                ALIGNER=aligner
            )
        )
        quit(1)

    pass


def extractUmappedReads(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# extract unmapped reads
#
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped_R1.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}_unmapped.fastq \\
        {PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Unmapped read extraction completed${{reset}}"
else
    logthis "Unmapped reads extracted to initial FASTQ, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            PIPELINE=pipeline,
        )
    )


def alignAndSort(script: TextIOWrapper, options: OptionsDict):
    processUnmapped = options["processUnmapped"]
    alignOnly = options["alignOnly"]

    filenames = getFileNames(options)
    trimmedFilenames = getTrimmedFileNames(options)

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filenames[0], trimmedFilenames[0], options)
    alignFASTQ(script, trimmedFilenames[0], options)
    sortAlignedAndMappedData(script, options)
    generateDepth(script, options)

    if processUnmapped == True:
        extractUmappedReads(script, options)

    if alignOnly == True:
        script.write(
            """

# align-only flag set
logthis "align-only set, exiting pipeline early"
exit
"""
        )


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
        choices=["bwa", "bwa-mem2", "hisat2"],
        help="Use 'bwa', 'bwa-mem2', or 'hisat2' as the aligner.",
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


def main(panel_choices: List[str], panel_choice_help: str):
    opts = defineArguments(panel_choices, panel_choice_help)
    options = fixupPathOptions(opts)
    verifyFileNames(options)
    verifyOptions(options)

    filenames = getFileNames(options)

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")

        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        updateDictionary(script, options, panel_choices)

        alignAndSort(script, options)
        runVariantPipeline(script, options)
        assignClade(script, options)

        # we'll wait here to make sure all the background stuff is done before we
        # run multiqc and cleanup
        script.write('logthis "${green}Done with front-end processing${reset}"\n')

        if options["runQc"] == True:
            doQualityControl(script, options, filenames)
            script.write('logthis "${green}Done with back-end processing${reset}"\n')

    system("chmod +x " + options["script"])
