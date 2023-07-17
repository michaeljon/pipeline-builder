#!/usr/bin/env python

import sys
import os
from io import TextIOWrapper
from argparse import ArgumentParser, Namespace
from os.path import exists, expandvars

from bio_types import *
from bio_methods import commonPipeline
from common import pipelineDriver
from single_end import runBwaAligner


def getFileNames(options: OptionsDict) -> FastqSet:
    fastqs = options["fastqs"]

    for fastq in fastqs:
        if exists(expandvars(fastq)) == False:
            print("Check your --fastqs parameters, some files were not found")
            quit(1)

    return fastqs


def verifyFileNames(options: OptionsDict):
    fastqs = getFileNames(options)

    for fastq in fastqs:
        if exists(expandvars(fastq)) == False:
            print("Check your --fastqs parameters, file {} was not found".format(fastq))
            quit(1)


def getTargetFilename(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return [
        "{PIPELINE}/{SAMPLE}.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    ]


def defineExtraArguments(parser: ArgumentParser):
    parser.add_argument(
        "--fastqs",
        nargs="+",
        action="store",
        metavar="FASTQS",
        dest="fastqs",
        default=[],
    )


def preprocess(script: TextIOWrapper, options: OptionsDict):
    fastqs = " ".join(getFileNames(options))
    fastq = getTargetFilename(options)

    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {FASTQ} ]]; then
    logthis "${{yellow}}Concatenating lane data${{reset}}"

    cat {FASTQS} >{FASTQ}

    logthis "${{yellow}}Lane data concatenated${{reset}}"
else
    logthis "Lane data already concatenated, ${{green}}skipping${{reset}}"
fi
""".format(
            FASTQS=fastqs,
            FASTQ=fastq[0],
        )
    )


def align(script: TextIOWrapper, options: OptionsDict, reference):
    fastq = getTargetFilename(options)

    runBwaAligner(script, fastq[0], options, reference)


def main():
    pipelineDriver(
        defineExtraArguments,
        verifyFileNames,
        getFileNames,
        lambda script, options, filenames, references: commonPipeline(
            script,
            options,
            "Illumina.SingleEnd",
            filenames,
            references,
            preprocess,
            align,
        ),
    )
