#!/usr/bin/env python

import argparse
import re
import csv
import gzip

from urllib import parse
from os.path import exists, expandvars

#
# list of primary regions in fasta
# grep -e '^>[^ ]* Homo sapiens chromosome [1-9XY]*, GRCh38.p14 Primary Assembly' GCF_000001405.40_GRCh38.p14_genomic.fna

# sanity-check regions in fasta
# awk '$1~/>/{a[$1]++}END{for (e in a) { if (a[e] > 1) { print e, a[e] }}}' GCF_000001405.40_GRCh38.p14_genomic.fna

# count of features per region
# awk '!($1~/#/){a[$1]++}END{for (e in a) { print e, a[e] }}' GCF_000001405.40_GRCh38.p14_genomic.gff  | sort -k2 -g -r


# RefSeq sequences with the prefixes NC_, NT_, NW_,NG_, NM_, NR_ or NP_
#   chromosome - NC_000023.11
#   genomic contigs or scaffolds - NT_010718.17, NW_003315950.2
#   gene/genomic region - NG_012232.1
#   coding transcript - NM_004006.2
#   non-coding transcript - NR_004430.2
#   protein - NP_003997.1

# index int8range layout

#                             bit position

#   6 6 6 6 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3
#   3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2
# [ p|  reserved   |  chromosome   |           sequence             ]

#   3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
#   1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0
# [                            locus                                ]
#

parser = argparse.ArgumentParser()
parser.add_argument(
    "-g",
    "--gff",
    required=True,
    action="store",
    metavar="GFF",
    dest="gff",
    help="Path to feature file",
)
parser.add_argument(
    "-f",
    "--fasta",
    required=True,
    action="store",
    metavar="FASTA",
    dest="fasta",
    help="Path to FASTA file",
)
parser.add_argument(
    "-r",
    "--ref",
    required=True,
    action="store",
    metavar="REFERENCE",
    dest="ref",
    help="Reference name to store",
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    action="store",
    metavar="TSV",
    dest="tsv",
    help="Path to output TSV",
)
opts = parser.parse_args()
options = vars(opts)

for opt in ["gff", "fasta", "tsv"]:
    options[opt] = expandvars(options[opt])

if exists(options["gff"]) == False:
    print("Unable to find GFF (--gff) at {PATH}".format(PATH=options["gff"]))
    quit(1)

if exists(options["fasta"]) == False:
    print("Unable to find FASTA (--fasta) at {PATH}".format(PATH=options["fasta"]))
    quit(1)

ISPRIMARY_MASK = 0x4000000000000000
RESERVED_MASK = 0xFF00000000000000
CHROMOSOME_MASK = 0x00FF000000000000
SEQUENCE_MASK = 0x0000FFFF00000000
LOCUS_MASK = 0x00000000FFFFFFFF

sequences = {}
sequenceTypes = {
    "NC": "chromosome",
    "NT": "genomic contig or scaffold",
    "NW": "genomic contig or scaffold",
    "NG": "gene/genomic region",
    "NM": "coding transcript",
    "NR": "non-coding transcript",
    "NP": "protein",
}


def processFeatures(writer):
    lineNumber = 0

    print("Processing features from GFF")

    gff = gzip.open(options["gff"], "rt") if options["gff"].endswith("gz") else open(options["gff"], "r")

    while line := gff.readline().rstrip():
        if line.startswith("#") == False:
            lineNumber += 1

            columns = line.split("\t")

            # if we've been asked to process a feature where there's no
            # matching sequence in the FASTA, we have no choice but to bail
            if not "sequenceTypeId" in sequences[columns[0]]:
                return

            sequenceTypeId = sequences[columns[0]]["sequenceTypeId"]

            sequence_id = columns[0]
            feature_type = columns[2]
            start_position = int(columns[3])
            stop_position = int(columns[4])

            info = columns[8].split(";")
            infos = {
                k: v for k, v in zip([k.split("=")[0] for k in info], [":".join(k.split("=")[1:]) for k in info])
            }

            id = infos["ID"] if "ID" in infos else None
            gene = infos["gene"] if "gene" in infos else None
            name = parse.unquote(infos["Name"]) if "Name" in infos else None
            description = parse.unquote(infos["description"]) if "description" in infos else None

            chromosome = sequences[columns[0]]["chromosome"]
            startKey = 0x0000000000000000

            if chromosome == "X":
                startKey |= 23 << 48
            elif chromosome == "Y":
                startKey |= 24 << 48
            # mitochondria
            elif chromosome == "M":
                startKey |= 25 << 48
            # unplaced
            elif chromosome == "Z":
                startKey |= 32 << 48
            # a numeric chromosome
            else:
                startKey |= int(chromosome) << 48

            startKey |= sequenceTypeId << 32
            if sequences[columns[0]]["isPrimary"] == True:
                startKey |= ISPRIMARY_MASK

            # endkey has identical high-order bits
            endKey = startKey

            # but the lower-order 32 bits are different
            startKey |= start_position
            endKey |= stop_position

            chromosome = sequences[columns[0]]["chromosomeName"]

            writer.writerow(
                [
                    options["ref"],
                    sequence_id,
                    feature_type,
                    start_position,
                    stop_position,
                    id,
                    name,
                    description,
                    gene,
                    startKey,
                    endKey,
                    None,
                    None,
                ]
            )

    gff.close()
    print(str(lineNumber) + " features processed from GFF")


def processSequence(line, sequenceId):
    rex = r">(?P<sequenceType>(NC|NT|NW|NG|NM|NR|NP))_(?P<sequence_id>[^ ]*) Homo sapiens (?P<flavor>(chromosome|unplaced.*|mitochondrion.*)) (?P<chromosome>[0-9XY]*)"
    sMatch = re.search(rex, line)

    if sMatch == None:
        return

    sequenceType = sMatch.group("sequenceType")
    sequence_id = sMatch.group("sequence_id")
    chromosome = sMatch.group("chromosome")
    flavor = sMatch.group("flavor").split(" ")[0]
    description = line[1:]

    sequence = {
        "sequenceTypeCode": sequenceType,
        "sequenceType": sequenceTypes[sequenceType],
        "sequenceTypeId": sequenceId,
        "sequence_id": sequenceType + "_" + sequence_id,
        "description": description,
        "isPrimary": False,
    }

    match flavor:
        case "unplaced":
            sequence["chromosomeName"] = "unplaced"
            sequence["chromosome"] = "Z"

        case "chromosome":
            sequence["chromosomeName"] = "chr" + chromosome
            sequence["chromosome"] = chromosome
            sequence["isPrimary"] = sequenceType == "NC"

        case "mitochondrion,":
            sequence["chromosomeName"] = "chrM"
            sequence["chromosome"] = "M"

        case _:
            sequence = {}

    sequences[sequenceType + "_" + sequence_id] = sequence


def processSequences():
    sequenceId = 1024

    print("Processing sequences from FASTA")

    fasta = gzip.open(options["fasta"], "rt") if options["fasta"].endswith("gz") else open(options["fasta"], "r")

    while line := fasta.readline().rstrip():
        if line.startswith(">"):
            sequenceId += 1
            processSequence(line, sequenceId)

    fasta.close()

    print(str(sequenceId - 1024) + " sequences processed from FASTA")


def main():
    processSequences()

    if options["tsv"] == "-":
        options["tsv"] = "/dev/stdout"

    with open(options["tsv"], "w") as tsv:
        writer = csv.writer(tsv, delimiter="\t")
        writer.writerow(
            [
                "reference",
                "sequence_id",
                "feature_type",
                "start_position",
                "stop_position",
                "id",
                "name",
                "description",
                "gene",
                "startKey",
                "endKey",
                "local_region",
                "global_region",
            ]
        )
        processFeatures(writer)


if __name__ == "__main__":
    main()
