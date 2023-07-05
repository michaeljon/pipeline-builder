from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars
from os import cpu_count, system

from bio_types import *
from common import *
from bio_methods import *

fallback_warning_shown = False


def getFileNames(options: OptionsDict) -> FastqSet:
    global fallback_warning_shown

    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    filenames = [
        options["r1"],
        options["r2"],
    ]

    if exists(expandvars(filenames[0])) == True and exists(expandvars(filenames[1])) == True:
        return filenames

    # assume we have the _001 pattern first
    filenames = [
        "{FASTQ_DIR}/{SAMPLE}_R1_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    ]

    if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
        if fallback_warning_shown == False:
            print("Falling back to shortened fastq file names")
            fallback_warning_shown = True

        # if that didn't work, try for the redacted names
        filenames = [
            "{FASTQ_DIR}/{SAMPLE}_R1.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
            "{FASTQ_DIR}/{SAMPLE}_R2.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        ]

        if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
            print(
                "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                    R1=filenames[0],
                    R2=filenames[1],
                )
            )
            print("Check your --sample and --fastq-dir parameters")
            quit(1)

    return filenames


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return [
        "{PIPELINE}/{SAMPLE}_R1.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    ]


def runIdentityPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
):
    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running identity preprocessor${{reset}}"

    ln -s {R1} {O1}
    ln -s {R2} {O2}

    logthis "${{yellow}}Identity preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
        )
    )


def runFastpPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    adapters = options["adapters"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    readLimit = int(options["read-limit"])

    #
    # fastp is limited to 8 threads, so we use them
    #

    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running FASTP preprocessor${{reset}}"

    fastp \\
        --report_title "fastp report for sample {SAMPLE}" \\
        --in1 {R1} \\
        --in2 {R2} \\
        --out1 {O1} \\
        --out2 {O2} \\
        {ADAPTERS} \\
        --verbose {LIMITREADS} \\
        --thread 8 \\
        -j {STATS}/{SAMPLE}-fastp.json \\
        -h {STATS}/{SAMPLE}-fastp.html

    logthis "${{yellow}}FASTP preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            ADAPTERS="--adapter_fasta " + adapters if adapters != "" else "--detect_adapter_for_pe",
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def runTrimmomaticPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    bin = options["bin"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]

    u1 = o1.replace(".trimmed", ".unpaired")
    u2 = o2.replace(".trimmed", ".unpaired")

    script.write(
        """
#
# run the trimmomatic preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running trimmomatic preprocessor${{reset}}"

    java -jar {BIN}/trimmomatic-0.39.jar \\
        PE \\
        -phred33 \\
        -threads {THREADS} \\
        {R1} \\
        {R2} \\
        {O1} {U1} \\
        {O2} {U2} \\
        ILLUMINACLIP:{BIN}/adapters/Ovation-PE.fa:2:30:5:2:True \\
        HEADCROP:15 \\
        LEADING:3 \\
        TRAILING:3 \\
        MINLEN:30 2> {STATS}/{SAMPLE}_trim_out.log

    logthis "${{yellow}}trimmomatic preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
            U1=u1,
            U2=u2,
            STATS=stats,
            SAMPLE=sample,
            BIN=bin,
            THREADS=threads,
        )
    )


def runTrimGalore(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    bin = options["bin"]
    sample = options["sample"]
    threads = options["cores"]
    pipeline = options["pipeline"]
    stats = options["stats"]

    #
    # trim_galore is limited to 8 threads, so we use them
    #

    script.write(
        """
#
# run the trim galore preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running trim galore preprocessor${{reset}}"

    trim_galore \\
        --cores 8 \\
        --paired \\
        --clip_R2 10 \\
        --basename {SAMPLE} \\
        --output_dir {PIPELINE} \\
        -a NNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \\
        -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA \\
        {R1} \\
        {R2} 2> {STATS}/{SAMPLE}_trim_out.log

    mv {PIPELINE}/{SAMPLE}_val_1.fq.gz {O1}
    mv {PIPELINE}/{SAMPLE}_val_2.fq.gz {O2}

    logthis "${{yellow}}trim galore preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
            BIN=bin,
            THREADS=threads,
            STATS=stats,
            SAMPLE=sample,
            PIPELINE=pipeline,
        )
    )


def preprocessFASTQ(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    preprocessor = options["preprocessor"]

    if preprocessor == "none":
        runIdentityPreprocessor(script, r1, r2, o1, o2)
    elif preprocessor == "fastp":
        runFastpPreprocessor(script, r1, r2, o1, o2, options)
    elif preprocessor == "trimmomatic":
        runTrimmomaticPreprocessor(script, r1, r2, o1, o2, options)
    elif preprocessor == "trimgalore":
        runTrimGalore(script, r1, r2, o1, o2, options)
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)

    pass


def runBwaAligner(
    script: TextIOWrapper,
    o1: str,
    o2: str,
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
        -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        -t {THREADS} \\
        -Y \\
        -M \\
        -v 1 \\
        {REFERENCE}/{ASSEMBLY}.fna \\
        {O1} \\
        {O2} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            ALIGNER=aligner,
            O1=o1,
            O2=o2,
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
    o1: str,
    o2: str,
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
            -1 {O1} \\
            -2 {O2} \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:ILLUMINA" \\
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
            O1=o1,
            O2=o2,
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
    o1: str,
    o2: str,
    options: OptionsDict,
):
    aligner = options["aligner"]

    if aligner == "bwa" or aligner == "bwa-mem2":
        runBwaAligner(script, o1, o2, options)
    elif aligner == "hisat2":
        runHisatAligner(script, o1, o2, options)
    else:
        print("Unexpected value {ALIGNER} given for the --aligner option".format(ALIGNER=aligner))
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
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped_R1.fastq || ! -f {PIPELINE}/{SAMPLE}_unmapped_R2.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}_unmapped_R1.fastq \\
        -2 {PIPELINE}/{SAMPLE}_unmapped_R2.fastq \\
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

    preprocessFASTQ(
        script,
        filenames[0],
        filenames[1],
        trimmedFilenames[0],
        trimmedFilenames[1],
        options,
    )
    alignFASTQ(script, trimmedFilenames[0], trimmedFilenames[1], options)
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
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>_R[12].fastq.gz",
    )

    parser.add_argument(
        "--r1",
        required=True,
        action="store",
        metavar="R1",
        dest="r1",
        help="Full path to the forward (R1) read FASTQ",
    )

    parser.add_argument(
        "--r2",
        required=True,
        action="store",
        metavar="R2",
        dest="r2",
        help="Full path to the reverse (R2) read FASTQ",
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
        choices=["trimmomatic", "fastp", "trimgalore", "none"],
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
    filenames = getFileNames(options)

    if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
        print(
            "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                R1=filenames[0],
                R2=filenames[1],
            )
        )
        print("Check your --sample and --work-dir parameters")
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
        script.write("which bash\n")

        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

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
