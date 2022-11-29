#!/usr/bin/env python

import argparse
import math
import gzip
import sys

from math import ceil
from os.path import exists, expandvars
from pprint import pprint
from sys import stderr
from typing import Dict, Tuple, List, Any

OptionsDict = Dict[str, Any]

VCF_HEADER = "##fileformat=VCFv4.2"


def get_options() -> OptionsDict:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--known-sites",
        required=True,
        type=str,
        action="store",
        dest="known-sites",
        help="Location of known-sites file",
    )

    parser.add_argument(
        "--features",
        required=True,
        type=str,
        action="store",
        dest="features",
        help="Location of feature (GFF) file",
    )

    # parser.add_argument(
    #     "--reference",
    #     required=True,
    #     type=str,
    #     action="store",
    #     dest="reference",
    #     help="Location of reference (FNA) file",
    # )

    parser.add_argument(
        "--output",
        required=True,
        type=str,
        action="store",
        dest="output",
        help="Full path to output directory for interval files",
    )

    parser.add_argument(
        "--skip-headers",
        action="store_true",
        dest="skipHeaders",
        default=False,
        help="Skip writing VCF headers (other than those required by VCF)",
    )

    parser.add_argument(
        "--process-set",
        action="store",
        dest="process-set",
        default="all",
        choices=["all", "primary", "secondary"],
        help="Generate files for specific contigs.",
    )

    parser.add_argument(
        "--segment",
        type=int,
        action="store",
        metavar="SEGMENT_SIZE",
        dest="segmentSize",
        default=50_000_000,
        help="Size of interval partition (ADVANCED)",
    )

    parser.add_argument(
        "--factor",
        type=float,
        action="store",
        metavar="FACTOR",
        dest="factor",
        default=0.25,
        help="Interval remainder buffer (between 0.10 and 0.50) (ADVANCED)",
    )

    parser.add_argument(
        "-l",
        "--limit",
        required=False,
        action="store",
        metavar="LIMIT",
        dest="limit",
        default=sys.maxsize,
        type=int,
        help="Limit number of input variants to process",
    )

    opts = parser.parse_args()
    return vars(opts)


def verifyOptions(options: OptionsDict):
    if exists(options["known-sites"]) == False:
        print("Unable to find your known-sites at {PATH}".format(PATH=options["known-sites"]))
        quit(1)

    if exists(options["features"]) == False:
        print("Unable to find your features at {PATH}".format(PATH=options["features"]))
        quit(1)

    # if exists(options["reference"]) == False:
    #     print("Unable to find your reference at {PATH}".format(PATH=options["reference"]))
    #     quit(1)

    if exists(options["output"]) == False:
        print("Unable to find your output directory {PATH}".format(PATH=options["output"]))
        quit(1)


def read_vcf_headers(options: OptionsDict):
    if options["skipHeaders"] == True:
        headers = [VCF_HEADER]
    else:
        headers = []

        file = (
            gzip.open(options["known-sites"], "rt")
            if options["known-sites"].endswith(".gz")
            else open(options["features"], "r")
        )

        while line := file.readline().strip():
            if line.startswith("##") == False:
                break
            else:
                headers.append(line)

        file.close()

    return headers


def write_vcf_header(options: OptionsDict, interval, headers: list[str]):
    full_path = options["output"] + "/" + interval["filename"] + ".vcf.gz"

    f = gzip.open(full_path, "wt")
    for header in headers:
        f.write(header + "\n")
    return f


def write_vcf_headers(options: OptionsDict, features, headers: list[str]):
    for feature, intervals in features.items():
        for interval in intervals:
            interval["filename"] = write_vcf_header(options, interval, headers)


def read_features(options: OptionsDict):
    features = {}
    segmentSize = options["segmentSize"]
    factor = options["factor"]
    lastBlockMax = math.floor(segmentSize * factor)

    file = (
        gzip.open(options["features"], "rt") if options["features"].endswith(".gz") else open(options["features"], "r")
    )

    while line := file.readline().rstrip():
        rule = line.split(" ")

        if line.startswith("##sequence-region") == False:
            continue

        sequence_id = rule[1]
        length = int(rule[3])

        # if this is a primary contig (NC_*) and we're not processing them, skip
        if not (options["process-set"] == "all" or options["process-set"] == "primary") and sequence_id.startswith(
            "NC_"
        ):
            continue

        # if this is a secondary contig and we're not processing them, skip
        if not (
            options["process-set"] == "all" or options["process-set"] == "secondary"
        ) and not sequence_id.startswith("NC_"):
            continue

        if not (sequence_id in features):
            features[sequence_id] = []

        remainder = length - segmentSize
        segments = ceil(length / segmentSize)
        segment = 0

        while remainder > lastBlockMax:
            lower = segment * segmentSize + 1
            upper = (segment + 1) * segmentSize
            filename = "{sequence_id}:{lower}-{upper}".format(
                sequence_id=sequence_id,
                lower=int(lower),
                upper=int(upper),
            )

            features[sequence_id].append(
                {
                    "sequence_id": sequence_id,
                    "lower": lower,
                    "upper": upper,
                    "filename": filename,
                    "file": None,
                }
            )

            segment += 1
            remainder -= segmentSize

        if remainder > 0:
            lower = (segments - 2) * segmentSize
            if lower % 10 == 0:
                lower += 1

            filename = "{sequence_id}:{lower}-{upper}".format(
                sequence_id=sequence_id,
                lower=lower,
                upper=length,
            )

            features[sequence_id].append(
                {
                    "sequence_id": sequence_id,
                    "lower": int(lower),
                    "upper": int(length),
                    "filename": filename,
                    "file": None,
                }
            )
        else:
            lower = (segments - 1) * segmentSize
            if lower % 10 == 0:
                lower += 1

            filename = "{sequence_id}:{lower}-{upper}".format(
                sequence_id=sequence_id,
                lower=lower,
                upper=length,
            )

            features[sequence_id].append(
                {
                    "sequence_id": sequence_id,
                    "lower": int(lower),
                    "upper": int(length),
                    "filename": filename,
                    "file": None,
                }
            )

    file.close()
    return features


def get_file_from_locus_in_features(sequence_id: str, position: int, features):
    if sequence_id not in features:
        return None

    chromosome = features[sequence_id]
    for interval in chromosome:
        if interval["lower"] <= position and position <= interval["upper"]:
            return interval["filename"]

    print("Unable to locate interval for {CHROM} at {POSITION}".format(CHROM=chromosome, POSITION=position))
    quit(1)


def process_known_sites(options: OptionsDict, features):
    cycle_indicators = ["|", "/", "-", "\\"]
    cycle = 0

    reads_remaining = options["limit"]

    ks = (
        gzip.open(options["known-sites"], "rt")
        if options["known-sites"].endswith(".gz")
        else open(options["features"], "r")
    )

    while line := ks.readline().rstrip():
        if reads_remaining <= 0:
            break

        if line.startswith("#") == False:
            columns = line.split("\t")
            sequence_id = columns[0]
            position = int(columns[1])

            # if this is a primary contig (NC_*) and we're not processing them, skip
            if not (options["process-set"] == "all" or options["process-set"] == "primary") and sequence_id.startswith(
                "NC_"
            ):
                continue

            # if this is a secondary contig and we're not processing them, skip
            if not (
                options["process-set"] == "all" or options["process-set"] == "secondary"
            ) and not sequence_id.startswith("NC_"):
                continue

            reads_remaining -= 1
            if reads_remaining % 20000 == 0:
                print("\b" + str(cycle_indicators[cycle % len(cycle_indicators)]), end="", flush=True)
                cycle += 1

            of = get_file_from_locus_in_features(sequence_id, position, features)
            if of != None:
                of.write(line + "\n")

    ks.close()
    print(" finished with " + str(options["limit"] - reads_remaining) + " processed", flush=True)


def main():
    options = get_options()
    verifyOptions(options)

    # read the chromosomes from the bam via stdin
    print("Reading features...", end="", flush=True)
    features = read_features(options)
    print(str(len(features.keys())) + " loaded", flush=True)

    # read the vcf headers from the known sites file
    print("Reading known site VCF headers...", end="", flush=True)
    headers = read_vcf_headers(options)
    print(str(len(headers)) + " loaded", flush=True)

    # write the headers and collect the output files
    print("Opening intervals and writing headers...", end="", flush=True)
    write_vcf_headers(options, features, headers)
    print(str(len(headers)) + " written to " + str(len(features.keys())) + " intervals", flush=True)

    # create intervals from known sites vcf
    print("Processing known sites... ", end="", flush=True)
    process_known_sites(options, features)

    # close all the open files
    print("Closing interval files", flush=True)
    for feature, intervals in features.items():
        for interval in intervals:
            interval["filename"].close()


if __name__ == "__main__":
    main()
