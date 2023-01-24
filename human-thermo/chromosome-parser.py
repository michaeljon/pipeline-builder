#!/usr/bin/env python

import json
import argparse
from os.path import exists, expandvars

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    required=True,
    action="store",
    metavar="GFF",
    dest="gff",
    help="Name of the 'feature' GFF",
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    action="store",
    metavar="JSON",
    dest="json",
    help="Name of the output 'chromosize size' file",
)
opts = parser.parse_args()
options = vars(opts)

for opt in ["gff", "json"]:
    options[opt] = expandvars(options[opt])

if exists(options["gff"]) == False:
    print("Unable to find GFF (--input) at {PATH}".format(PATH=options["gff"]))
    quit(1)

accessions = []


def processFeature(line):
    columns = line.split("\t")

    if columns[1] != "RefSeq" or columns[2] != "region" or not "genome=chromosome" in columns[8]:
        return

    accession = columns[0]
    length = int(columns[4])
    infos = columns[8].split(";")
    info2 = {k: v for k, v in zip([k.split("=")[0] for k in infos], [":".join(k.split("=")[1:]) for k in infos])}

    accessions.append({"chromosome": info2["Name"], "length": length, "accession": accession})


gff = options["gff"]

print("Processing features from GFF")
with open(gff, "r") as gff:
    while line := gff.readline().rstrip():
        if line.startswith("#"):
            continue

        processFeature(line)

print("GFF processed")
with open(options["json"], "w") as sizes:
    json.dump(accessions, sizes, indent=2, sort_keys=False)
