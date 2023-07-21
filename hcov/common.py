from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars
from datetime import datetime
from os import cpu_count, system, path, mkdir
import sys
import json

from bio_types import *


def writeVersions(script: TextIOWrapper):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...\n")
    script.write("#\n")


def writeEnvironment(script: TextIOWrapper, options: OptionsDict):
    noColor = options["noColor"]

    if noColor == True:
        script.write(
            """
# color shortcuts disabled
red=$(echo -n "")
green=$(echo -n "")
yellow=$(echo -n "")
reset=$(echo -n "")
"""
        )
    else:
        script.write(
            """
# color shortcuts enabled
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
reset=$(tput sgr0)
"""
        )

    script.write("#\n")
    script.write(
        """
CURRENT_REFERENCE='*';
function logthis() {{
  NOW=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "[${{NOW}}] (${{CURRENT_REFERENCE}}) ${{1}}"
}}

# deal with really large fastq files
ulimit -n 8192

# perl stuff
export PATH={WORKING}/perl5/bin:{WORKING}/bin/gatk-4.3.0.0:$PATH
export PERL5LIB={WORKING}/perl5/lib/perl5:$PERL5LIB
export PERL_LOCAL_LIB_ROOT={WORKING}/perl5:$PERL_LOCAL_LIB_ROOT

# shared library stuff (linux)
export LD_LIBRARY_PATH={WORKING}/lib:{WORKING}/bin:/usr/lib64:/usr/local/lib/:$LD_LIBRARY_PATH

# shared library stuff (darwin)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:{WORKING}/lib:{WORKING}/bin

# bcftools
export BCFTOOLS_PLUGINS={WORKING}/libexec/bcftools

# handy path
export PATH={WORKING}/bin/FastQC:{WORKING}/bin:$PATH\n""".format(
            WORKING=options["working"]
        )
    )
    script.write("\n")


def fixupPathOptions(opts: Namespace) -> OptionsDict:
    options = vars(opts)

    if options["sample"] == None:
        print("--sample is a required option")
        quit(1)

    if options["working"] == None:
        options["working"] = "$HOME"

    if options["bin"] == None:
        options["bin"] = "{WORKING}/bin".format(WORKING=options["working"])
    if options["lib"] == None:
        options["lib"] = "{WORKING}/lib".format(WORKING=options["working"])
    if options["reference"] == None:
        options["reference"] = "{WORKING}/reference".format(WORKING=options["working"])
    if options["temp"] == None:
        options["temp"] = "{WORKING}/temp".format(WORKING=options["working"])

    if options["pipeline"] == None:
        options["pipeline"] = "{WORKING}/pipeline/{SAMPLE}".format(
            WORKING=options["working"],
            SAMPLE=options["sample"],
        )
    if options["stats"] == None:
        options["stats"] = "{WORKING}".format(
            WORKING=options["pipeline"],
        )

    if options["script"] == None:
        options["script"] = "{OUTPUT}/{SAMPLE}_runner".format(
            OUTPUT=options["pipeline"],
            SAMPLE=options["sample"],
        )

    for opt in [
        "working",
        "reference",
        "pipeline",
        "stats",
        "temp",
        "bin",
        "script",
        "adapters",
    ]:
        options[opt] = expandvars(options[opt])

    return options


def defineArguments(defineExtraArguments) -> Namespace:
    parser = ArgumentParser()

    parser.add_argument(
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="Short name of sample, all files will use this as an alias",
    )

    if defineExtraArguments != None:
        defineExtraArguments(parser)

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
        "--reference-rules",
        action="store",
        metavar="REFERENCE_RULES",
        dest="referenceRules",
        help="Full path to the references.json file holding accessing and assembly rules",
        default=path.join(path.dirname(sys.argv[0]), "..", "references.json"),
    )

    parser.add_argument(
        "--references",
        nargs="+",
        required=True,
        action="store",
        metavar="REFERENCES",
        dest="references",
    )

    # then this becomes the root of the reference dir
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


def verifyOptions(options: OptionsDict):
    if exists(options["bin"]) == False:
        print("Unable to find your --bin-dir directory at {PATH}".format(PATH=options["bin"]))
        quit(1)

    if exists(options["lib"]) == False:
        print("Unable to find your --lib-dir directory at {PATH}".format(PATH=options["lib"]))
        quit(1)

    if exists(options["working"]) == False:
        print("Unable to find your --work-dir directory at {PATH}".format(PATH=options["working"]))
        quit(1)

    if exists(options["temp"]) == False:
        print("Unable to find your --temp-dir directory at {PATH}".format(PATH=options["temp"]))
        quit(1)

    if exists(options["reference"]) == False:
        print("Unable to find your --reference-dir directory at {PATH}".format(PATH=options["reference"]))
        quit(1)

    if exists(options["pipeline"]) == False:
        mkdir(options["pipeline"])

    if exists(options["stats"]) == False:
        mkdir(options["stats"])


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: List[str]):
    script.write("#\n")
    script.write("# generated at {TIME}\n".format(TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    script.write("#\n")
    script.write("# Parameters\n")
    for opt in options.keys():
        script.write("#   {OPTION} = {VALUE}\n".format(OPTION=opt, VALUE=options[opt]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Input files\n")
    for i, f in enumerate(filenames):
        script.write("#   R{R} = {F}\n".format(F=f, R=i))
    script.write("#\n")


def pipelineDriver(
    defineExtraArguments,
    verifyFileNames,
    getFileNames,
    specificPipelineRunner,
):
    opts = defineArguments(defineExtraArguments)
    options = fixupPathOptions(opts)
    verifyFileNames(options)
    verifyOptions(options)

    references = {}
    with open(options["referenceRules"]) as f:
        references = json.load(f)

    for ref in references.keys():
        reference = references[ref]

        # stuff the key into the object so we have it
        reference["common"] = ref

        # and make sure we have a path to the file
        reference["assembly"] = path.join(options["reference"], reference["assembly"])
        if exists(reference["assembly"]) == False:
            print("Missing reference FNA at {PATH} for {REF}".format(REF=ref, PATH=reference["assembly"]))
            quit(1)

    filenames = getFileNames(options)

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")

        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        specificPipelineRunner(script, options, filenames, references)

    system("chmod +x " + options["script"])
