#!/usr/bin/env python

from hcov_oc43_common import (
    OptionsDict,
    updateDictionary,
    sortAlignedAndMappedData,
    fixupPathOptions,
    verifyOptions,
    writeVersions,
    writeEnvironment,
    runAlignmentQC,
    runPipeline,
    doVariantQC,
    runMultiQC,
)

from io import TextIOWrapper
import argparse

from argparse import Namespace
from typing import Tuple
from datetime import datetime
from os.path import exists, expandvars
from os import cpu_count, system

FastqSet = Tuple[str, str]


def getFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_R1_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    )


def runIdentityPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    cp {R1} {O1}
    cp {R2} {O2}
else
    echo "Preprocessor already run, ${{green}}skipping${{reset}}"
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

    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
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
else
    echo "Preprocessor already run, ${{green}}skipping${{reset}}"
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
    timeout = options["watchdog"]
    bin = options["bin"]

    script.write(
        """
#
# run the trimmomatic preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
    java -jar ~/bin/trimmomatic-0.39.jar PE \\
        {R1} \\
        {R2} \\
        {O1} /dev/null \\
        {O2} /dev/null \\
        ILLUMINACLIP:{BIN}/adapters/NexteraPE-PE.fa:2:30:10 \\
        LEADING:5 \\
        TRAILING:5 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:30
else
    echo "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1, R2=r2, O1=o1, O2=o2, TIMEOUT=timeout, BIN=bin
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
        runIdentityPreprocessor(script, r1, r2, o1, o2, options)
    elif preprocessor == "fastp":
        runFastpPreprocessor(script, r1, r2, o1, o2, options)
    elif preprocessor == "trimmomatic":
        runTrimmomaticPreprocessor(script, r1, r2, o1, o2, options)
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
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    nonRepeatable = options["non-repeatable"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    bwa-mem2 mem -t {THREADS} \\
        -Y -M {DASHK} \\
        -v 1 \\
        -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        {REFERENCE}/hcov-oc43.fasta \\
        {O1} \\
        {O2} | 
    samtools view -Sb - >{PIPELINE}/{SAMPLE}.aligned.bam
else
    echo "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            DASHK="" if nonRepeatable == True else "-K " + str((10_000_000 * int(threads))),
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
    hisat2 \\
            -x {REFERENCE}/hcov-oc43 \\
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
    samtools view -Sb - >{PIPELINE}/{SAMPLE}.aligned.bam
else
    echo "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            O1=o1,
            O2=o2,
            REFERENCE=reference,
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

    if aligner == "bwa":
        runBwaAligner(script, o1, o2, options)
    elif aligner == "hisat2":
        runHisatAligner(script, o1, o2, options)
    else:
        print("Unexpected value {ALIGNER} given for the --aligner option".format(ALIGNER=aligner))
        quit(1)

    pass


def alignAndSort(script: TextIOWrapper, options: OptionsDict, output: str):
    filenames = getFileNames(options)
    trimmedFilenames = getTrimmedFileNames(options)

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filenames[0], filenames[1], trimmedFilenames[0], trimmedFilenames[1], options)
    alignFASTQ(script, trimmedFilenames[0], trimmedFilenames[1], options)

    sortAlignedAndMappedData(script, options, output)


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastqSet):
    script.write("#\n")
    script.write("# generated at {TIME}\n".format(TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    script.write("#\n")
    script.write("# Parameters\n")
    for opt in options.keys():
        script.write("#   {OPTION} = {VALUE}\n".format(OPTION=opt, VALUE=options[opt]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Input files\n")
    script.write("#   R1 = {F}\n".format(F=filenames[0]))
    script.write("#   R2 = {F}\n".format(F=filenames[1]))
    script.write("#\n")


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.set_defaults(doQC=False, cleanIntermediateFiles=True)
    parser.add_argument(
        "-s",
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>_R[12].fastq.gz",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        required=True,
        action="store",
        metavar="WORKING_DIR",
        dest="working",
        help="Working directory, e.g. base for $WORKING/pipeline, $WORKING/stats",
    )
    parser.add_argument(
        "-q",
        "--skip-qc",
        action="store_false",
        dest="doQC",
        default=True,
        help="Skip QC process on input and output files",
    )
    parser.add_argument(
        "-P",
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["trimmomatic", "fastp", "none"],
        help="Optionally run a FASTQ preprocessor",
    )

    parser.add_argument(
        "-a",
        "--align-only",
        action="store_true",
        dest="alignOnly",
        default=False,
        help="Only run alignment and sorting processes",
    )

    parser.add_argument(
        "-V",
        "--caller",
        action="store",
        dest="caller",
        default="bcftools",
        choices=["bcftools", "gatk", "lofreq"],
        help="Use `bcftools`, `gatk`, or `lofreq` as the variant caller",
    )

    parser.add_argument(
        "-I",
        "--alternate-consensus",
        action="store_true",
        dest="alternateConsensus",
        default=False,
        help="Use alternate consensus generation process",
    )

    parser.add_argument(
        "-A",
        "--aligner",
        action="store",
        dest="aligner",
        default="bwa",
        choices=["bwa", "hisat2"],
        help="Use 'bwa' or 'hisat2' as the aligner.",
    )

    parser.add_argument(
        "-S",
        "--sorter",
        action="store",
        dest="sorter",
        default="biobambam",
        choices=["biobambam", "samtools"],
        help="Use 'biobambam' or 'samtools' as the sorter.",
    )

    parser.add_argument(
        "-r",
        "--reference-dir",
        action="store",
        metavar="REFERENCE_DIR",
        dest="reference",
        help="Location of the reference genome files",
    )
    parser.add_argument(
        "-p",
        "--pipeline-dir",
        action="store",
        metavar="PIPELINE_DIR",
        dest="pipeline",
        help="Output directory for processing",
    )
    parser.add_argument(
        "-F",
        "--fastq-dir",
        action="store",
        metavar="FASTQ_DIR",
        dest="fastq_dir",
        help="Location of R1 and R2 files",
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        action="store",
        metavar="TEMP_DIR",
        dest="temp",
        help="Temporary storage for alignment",
    )
    parser.add_argument(
        "-o",
        "--stats-dir",
        action="store",
        metavar="STATS_DIR",
        dest="stats",
        help="Destination for statistics and QC files",
    )
    parser.add_argument(
        "-b",
        "--bin-dir",
        action="store",
        metavar="BIN_DIR",
        dest="bin",
        default="$HOME/bin",
        help="Install location of all tooling",
    )

    parser.add_argument(
        "-D",
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
        "-z",
        "--cores",
        action="store",
        dest="cores",
        metavar="CPU_COUNT",
        default=cpu_count(),
        help="Specify the number of available CPU",
    )

    parser.add_argument(
        "-N",
        "--non-repeatable",
        action="store_true",
        dest="non-repeatable",
        default=False,
        help="Force the -K option to BWA, will use 10,000,000bp per thread",
    )

    parser.add_argument(
        "--read-limit",
        action="store",
        dest="read-limit",
        metavar="READ_LIMIT",
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    parser.add_argument(
        "-d",
        "--watchdog",
        action="store",
        metavar="WATCHDOG_TIMEOUT",
        dest="watchdog",
        default=150,
        help="Specify a watchdog timeout for the alignment process. Value is in minutes",
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


def main():
    opts = defineArguments()
    options = fixupPathOptions(opts)
    verifyFileNames(options)
    verifyOptions(options)

    filenames = getFileNames(options)
    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])
    prefix = "{PIPELINE}/{SAMPLE}".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])
    aligned = "{PIPELINE}/{SAMPLE}.aligned.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        script.write("\n")
        # script.write("set -e\n")
        # script.write("set -o pipefail\n")
        # script.write("\n")

        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write("touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

        updateDictionary(script, options)
        filenames = getFileNames(options)

        alignAndSort(
            script,
            options,
            sorted,
        )

        if options["doQC"]:
            runAlignmentQC(script, options, sorted, aligned)

        runPipeline(script, options, prefix)

        if options["doQC"]:
            script.write(
                """
echo ${yellow}Waiting for variant calling pipeline to complete before starting variant qc and multiqc${reset}
wait
echo ${green}Pipeline processed${reset}
                """
            )

            doVariantQC(script, options)
            runMultiQC(script, options)

        script.write(
            """
\necho ${{yellow}}Waiting for any outstanding processes to complete, this might return immediately, it might not.${{reset}}
wait
wait
\necho -e "${{green}}Done processing${{reset}} {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"\n""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                PIPELINE=options["pipeline"],
            )
        )

        script.write("\n")
        script.write("touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

    system("chmod +x " + options["script"])


if __name__ == "__main__":
    main()
