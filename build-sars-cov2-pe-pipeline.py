#!/usr/bin/env python

from sars_common import (
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
from typing import  Tuple
from datetime import datetime
from os.path import exists, expandvars
from os import cpu_count, system

FastqSet = Tuple[str, str]


def getFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_R1.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )


def alignWithFastpAndBWA(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    options: OptionsDict,
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    nonRepeatable = options["non-repeatable"]
    readLimit = int(options["read-limit"])

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --report_title "fastp report for sample {SAMPLE}" \\
            --in1 {R1} \\
            --in2 {R2} \\
            --detect_adapter_for_pe \\
            --merge \\
            --verbose {LIMITREADS} \\
            --stdout \\
            --thread 8 \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html | \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/covid_reference.fasta \\
            - | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            DASHK=""
            if nonRepeatable == True
            else "-K " + str((10_000_000 * int(threads))),
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )

def alignWithFastpAndHisat(
    script: TextIOWrapper, 
    r1: str,
    r2: str,
    options: OptionsDict
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    readLimit = int(options["read-limit"])

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --report_title "fastp report for sample {SAMPLE}" \\
            --in1 {R1} \\
            --in2 {R2} \\
            --merge \\
            --verbose {LIMITREADS} \\
            --stdout \\
            --thread 8 \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html |
        hisat2 \\
            -x {REFERENCE}/covid_reference \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:ILLUMINA" \\
            --rg "PU:unspecified" \\
            --rg "LB:{SAMPLE}" \\
            --rg "SM:{SAMPLE}" \\
            -U- \\
            --no-spliced-alignment \\
            --no-unal \\
            --threads {THREADS} | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )

def alignWithoutFastpUsingBWA(
    script: TextIOWrapper, r1: str, r2: str, options: OptionsDict
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    nonRepeatable = options["non-repeatable"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/covid_reference.fasta \\
            {R1} \\
            {R2} \\
            | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            DASHK=""
            if nonRepeatable == True
            else "-K " + str((10_000_000 * int(threads))),
            BIN=bin,
        )
    )

def alignWithoutFastpAndHisat(
    script: TextIOWrapper, 
    r1: str,
    r2: str,
    options: OptionsDict
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        hisat2 \\
            -x {REFERENCE}/covid_reference \\
            -1 {R1} \\
            -2 {R2} \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:ILLUMINA" \\
            --rg "PU:unspecified" \\
            --rg "LB:{SAMPLE}" \\
            --rg "SM:{SAMPLE}" \\
            --no-spliced-alignment \\
            --no-unal \\
            --threads {THREADS} | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            BIN=bin,
        )
    )

def alignAndSort(
    script: TextIOWrapper, options: OptionsDict, output: str
):
    skipFastp = options["skipPreprocess"]
    aligner = options["aligner"]
    filenames = getFileNames(options)

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    if skipFastp == False:
        if aligner == "bwa":
            alignWithFastpAndBWA(script, filenames[0], filenames[1], options)
        else:
            alignWithFastpAndHisat(script, filenames[0], filenames[1], options)
    else:
        if aligner == "bwa":
            alignWithoutFastpUsingBWA(script, filenames[0], filenames[1], options)
        else:
            alignWithoutFastpAndHisat(script, filenames[0], filenames[1], options)

    sortAlignedAndMappedData(script, options, output)


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastqSet):
    script.write("#\n")
    script.write(
        "# generated at {TIME}\n".format(
            TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
    )
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
        "--skip-fastp",
        action="store_true",
        dest="skipPreprocess",
        default=False,
        help="Skip running fastp on input file(s)",
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
        help="Use `bcftools` or `gatk` as the variant caller",
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
        help="Use 'bwa' or 'hisat2' as the aligner.",
    )

    parser.add_argument(
        "-S",
        "--sorter",
        action="store",
        dest="sorter",
        default="biobambam",
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

    if (
        exists(expandvars(filenames[0])) == False
        or exists(expandvars(filenames[1])) == False
    ):
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
    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )
    prefix = "{PIPELINE}/{SAMPLE}".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )
    aligned = "{PIPELINE}/{SAMPLE}.aligned.sam.gz".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write(
            "touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"])
        )
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
\necho -e "${{green}}Done processing${{reset}} {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"\n""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                PIPELINE=options["pipeline"],
            )
        )

        script.write("\n")
        script.write(
            "touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"])
        )
        script.write("\n")

    system("chmod +x " + options["script"])


if __name__ == "__main__":
    main()
