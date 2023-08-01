from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars


from bio_types import *
from bio_methods import *
from common import pipelineDriver


def getFileNames(options: OptionsDict) -> FastqSet:
    fastqs = [
        options["r1"],
        options["r2"],
    ]

    if exists(expandvars(fastqs[0])) == True and exists(expandvars(fastqs[1])) == True:
        return fastqs

    print("Check your --r1 and --r2 parameters, files were not found")
    quit(1)


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
    sample = options["sample"]
    adapters = options["adapters"]
    stats = options["stats"]
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
            ADAPTERS="--adapter_fasta " + adapters if adapters != "" else "--detect_adapter_for_pe",
            SAMPLE=sample,
            STATS=stats,
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
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)

    pass


def runBwaAligner(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    options: OptionsDict,
    reference,
):
    aligner = options["aligner"]
    referenceAssembly = reference["assembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    bin = options["bin"]
    instrument = options["_instrument"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {ALIGNED} ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    {ALIGNER} mem \\
        -R "@RG\\tID:{SAMPLE}\\tPL:{INSTRUMENT}\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        -t {THREADS} \\
        -v 1 \\
        {REFERENCE_ASSEMBLY} \\
        {R1} \\
        {R2} | 
    samtools view -Sb -@ 4 - >{ALIGNED}

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{ALIGNED}, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            ALIGNER=aligner,
            R1=r1,
            R2=r2,
            INSTRUMENT=instrument,
            REFERENCE_ASSEMBLY=referenceAssembly,
            ORGANISM=reference["common"],
            SAMPLE=sample,
            ALIGNED=buildBamFilePath(options, reference, "aligned"),
            THREADS=threads,
            PIPELINE=pipeline,
            BIN=bin,
        )
    )


def preprocess(script: TextIOWrapper, options: OptionsDict):
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


def align(script: TextIOWrapper, options: OptionsDict, reference):
    trimmedFilenames = getTrimmedFileNames(options)

    runBwaAligner(script, trimmedFilenames[0], trimmedFilenames[1], options, reference)


def verifyFileNames(options: OptionsDict):
    fastqs = getFileNames(options)

    if exists(expandvars(fastqs[0])) == False or exists(expandvars(fastqs[1])) == False:
        print(
            "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                R1=fastqs[0],
                R2=fastqs[1],
            )
        )
        print("Check your --r1 and --r2 parameters")
        quit(1)


def defineExtraArguments(parser: ArgumentParser):
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
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["trimmomatic", "fastp", "none"],
        help="Optionally run a FASTQ preprocessor",
    )


def main():
    pipelineDriver(
        defineExtraArguments,
        verifyFileNames,
        getFileNames,
        lambda script, options, filenames, references: commonPipeline(
            script,
            options,
            "Illumina.PairedEnd",
            filenames,
            references,
            preprocess,
            align,
        ),
    )
