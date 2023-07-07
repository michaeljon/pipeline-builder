from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars
from os import cpu_count

from bio_types import *
from bio_methods import *
from common import pipelineDriver


def getFileNames(options: OptionsDict) -> FastqSet:
    fastq = options["fastq"]

    if exists(expandvars(fastq)) == True:
        return [fastq]

    print("Check your --fastq parameter, files were not found")
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


def verifyFileNames(options: OptionsDict):
    fastq = options["fastq"]

    if exists(expandvars(fastq)) == True:
        return fastq

    print("Check your --fastq parameter")
    quit(1)


def defineExtraArguments(parser: ArgumentParser):
    parser.add_argument(
        "--fastq",
        required=True,
        action="store",
        metavar="FASTQ",
        dest="fastq",
        help="Full path to the FASTQ",
    )

    parser.add_argument(
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["cutadapt", "none"],
        help="Optionally run a FASTQ preprocessor",
    )


def main(panel_choices: List[str], panel_choice_help: str):
    pipelineDriver(
        panel_choices,
        panel_choice_help,
        defineExtraArguments,
        verifyFileNames,
        getFileNames,
        lambda s, o, f: commonPipeline(s, o, f, preprocessAndAlign),
    )
