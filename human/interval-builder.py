#!/usr/bin/env python

import argparse
import math

from math import ceil
from sys import stderr

parser = argparse.ArgumentParser()

parser.add_argument(
    "--process",
    type=str,
    required=True,
    action="store",
    metavar="PROCESS",
    dest="process",
    choices=["intervalList", "mergeList"],
    help="Generate an interval list or a VCF merge list from samtools --header-only input",
)

parser.add_argument(
    "--primaryOnly",
    action="store_true",
    dest="primaryOnly",
    default=False,
    help="Generate only the 'primary' sequences NC_*",
)

parser.add_argument(
    "--root",
    type=str,
    required=False,
    action="store",
    metavar="ROOT",
    dest="root",
    help="Root filename segment for VCF cleanup",
)

parser.add_argument(
    "--input",
    type=str,
    action="store",
    metavar="INPUT",
    dest="input",
    default="/dev/stdin",
    help="Input file holding list of @SQ headers",
)

parser.add_argument(
    "--output",
    type=str,
    action="store",
    metavar="OUTPUT",
    dest="output",
    default="/dev/stdout",
    help="Output interval file",
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

opts = parser.parse_args()
options = vars(opts)

if options["process"] == "mergeList" and options["root"] == None:
    stderr.write("When processing a merge list the VCF root component must be provided")
    exit(99)

intervals = []

segmentSize = options["segmentSize"]
factor = options["factor"]
lastBlockMax = math.floor(segmentSize * factor)

with open("/dev/stdin", "r") as file:
    while line := file.readline().rstrip():
        rule = line.split("\t")

        accession = rule[1].replace("SN:", "")
        length = int(rule[2].replace("LN:", ""))

        if options["primaryOnly"] == True and not accession.startswith("NC_"):
            continue

        remainder = length - segmentSize
        segments = ceil(length / segmentSize)
        segment = 0

        while remainder > lastBlockMax:
            lower = segment * segmentSize + 1
            upper = (segment + 1) * segmentSize
            intervals.append("{accession}:{lower}-{upper}".format(accession=accession, lower=lower, upper=upper))

            segment += 1
            remainder -= segmentSize

        if remainder > 0:
            lower = (segments - 2) * segmentSize
            if lower % 10 == 0:
                lower += 1

            intervals.append(
                "{accession}:{lower}-{upper}".format(
                    accession=accession,
                    lower=lower,
                    upper=length,
                )
            )
        else:
            lower = (segments - 1) * segmentSize
            if lower % 10 == 0:
                lower += 1

            intervals.append(
                "{accession}:{lower}-{upper}".format(
                    accession=accession,
                    lower=lower,
                    upper=length,
                )
            )

sorted_intervals = sorted(intervals, key=lambda k: (k.split(':')[0], int(k.split(':')[1].split('-')[0])))

with open("/dev/stdout", "w") as f:
    if options["process"] == "intervalList":
        f.write("interval\troot\n")

        for i in [(interval, interval.replace(":", "_").replace("-", "_")) for interval in sorted_intervals]:
            f.write("{INTERVAL}\t{ROOT}\n".format(INTERVAL=i[0], ROOT=i[1]))
    elif options["process"] == "mergeList": 
        for i in [interval.replace(":", "_").replace("-", "_") for interval in sorted_intervals]:
            f.write("{ROOT}.{INTERVAL}.vcf\n".format(ROOT=options["root"], INTERVAL=i))
    else:
        stderr.write("Invalid --process option??!!")
        exit(99)
