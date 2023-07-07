from io import TextIOWrapper
from argparse import Namespace
from os.path import exists, expandvars
from datetime import datetime
from os import cpu_count, system

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
function logthis() {{
  NOW=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "[${{NOW}}] ${{1}}"
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

    if options["working"] == None:
        options["working"] = "$HOME"
    if options["bin"] == None:
        options["bin"] = "{WORKING}/bin".format(WORKING=options["working"])
    if options["lib"] == None:
        options["lib"] = "{WORKING}/lib".format(WORKING=options["working"])
    if options["reference"] == None:
        options["reference"] = "{WORKING}/reference".format(WORKING=options["working"])
    if options["pipeline"] == None:
        options["pipeline"] = "{WORKING}/pipeline".format(WORKING=options["working"])
    if options["fastq_dir"] == None:
        options["fastq_dir"] = "{WORKING}/fastq".format(WORKING=options["working"])
    if options["stats"] == None:
        options["stats"] = "{WORKING}/stats".format(WORKING=options["working"])
    if options["temp"] == None:
        options["temp"] = "{WORKING}/temp".format(WORKING=options["working"])

    if options["sample"] == None:
        print("--sample is a required option")
        quit(1)

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
        print("Unable to find your --pipeline-dir directory at {PATH}".format(PATH=options["pipeline"]))
        quit(1)

    if exists(options["stats"]) == False:
        print("Unable to find your --stats-dir directory at {PATH}".format(PATH=options["stats"]))
        quit(1)

    ## this is a bit of a hack, but necessary since the tools and databases
    ## to do clade assignment and variant annotation don't exist right now
    referenceName = options["referenceName"]
    options["__canAssignClades"] = referenceName == "MN908947.3"
    options["__canAnnotateVariants"] = referenceName == "MN908947.3"


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
    panel_choices: List[str],
    panel_choice_help: str,
    defineArguments,
    verifyFileNames,
    getFileNames,
    specificPipelineRunner,
):
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

        specificPipelineRunner(script, options, filenames)

    system("chmod +x " + options["script"])
