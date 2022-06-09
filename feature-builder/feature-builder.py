from curses.ascii import isupper
import re
import json
from urllib import parse

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

FASTA = "/Users/michaeljon/Downloads/grch38p14 genome_assemblies_genome_gff/ncbi-genomes-2022-06-06/GCF_000001405.40_GRCh38.p14_genomic.fna"
GFF = "/Users/michaeljon/Downloads/grch38p14 genome_assemblies_genome_gff/ncbi-genomes-2022-06-06/GCF_000001405.40_GRCh38.p14_genomic.gff"

REFERENCE = "GRCh38.p14"
ACCESSION = "NCBI_Assembly:GCF_000001405.40"

ISPRIMARY_MASK = 0x8000000000000000
RESERVED_MASK = 0x7F00000000000000
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


def processFeature(line, lineNumber):
    feature = {}
    columns = line.split("\t")

    # if we've been asked to process a feature where there's no
    # matching sequence in the FASTA, we have no choice but to bail
    if ("sequenceTypeId" in sequences[columns[0]]) == False:
        return

    feature["featureLine"] = lineNumber
    feature["sequenceTypeId"] = sequences[columns[0]]["sequenceTypeId"]

    feature["reference"] = REFERENCE
    feature["assembly"] = ACCESSION
    feature["sequenceName"] = columns[0]

    feature["featureType"] = columns[2]
    feature["startLocus"] = int(columns[3])
    feature["endLocus"] = int(columns[4])
    feature["info"] = {}

    infos = columns[8].split(";")
    if infos != None and len(infos) > 0:
        for i in infos:
            v = i.split("=")
            if i.startswith("ID="):
                feature["id"] = v[1]
            elif i.startswith("gene="):
                feature["Gene"] = v[1]
            elif i.startswith("Dbxref="):
                identities = v[1].split(",")
                feature["alternateIdentities"] = {
                    k: v
                    for k, v in zip(
                        [k.split(":")[0] for k in identities], [":".join(k.split(":")[1:]) for k in identities]
                    )
                }
            elif isupper(i[0][0]):
                feature[v[0]] = v[1]
            else:
                feature["info"][v[0].lower()] = parse.unquote(v[1])

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

    startKey |= feature["sequenceTypeId"] << 32
    if sequences[columns[0]]["isPrimary"] == True:
        startKey |= ISPRIMARY_MASK

    # endkey has identical high-order bits
    endKey = startKey

    # but the lower-order 32 bits are different
    startKey |= feature["startLocus"]
    endKey |= feature["endLocus"]

    feature["chromosome"] = sequences[columns[0]]["chromosomeName"]

    feature["key"] = {"start": startKey, "end": endKey}
    feature["formattedKey"] = {"start": "0x" + format(startKey, "016x"), "end": "0x" + format(endKey, "016x")}

    print(json.dumps(feature, indent=2, sort_keys=False))
    # print(feature)


def processFeatures():
    lineNumber = 0

    print("Processing features from GFF")

    with open(GFF, "r") as gff:
        while line := gff.readline().rstrip():
            if line.startswith("#") == False:
                lineNumber += 1
                processFeature(line, lineNumber)

    print(str(lineNumber) + " features processed from GFF")


def processSequence(line, sequenceId):
    rex = r">(?P<sequenceType>(NC|NT|NW|NG|NM|NR|NP))_(?P<sequenceName>[^ ]*) Homo sapiens (?P<flavor>(chromosome|unplaced.*|mitochondrion.*)) (?P<chromosome>[0-9XY]*)"
    sMatch = re.search(rex, line)

    if sMatch == None:
        return

    sequenceType = sMatch.group("sequenceType")
    sequenceName = sMatch.group("sequenceName")
    chromosome = sMatch.group("chromosome")
    flavor = sMatch.group("flavor").split(" ")[0]
    description = line[1:]

    sequence = {
        "sequenceTypeCode": sequenceType,
        "sequenceType": sequenceTypes[sequenceType],
        "sequenceTypeId": sequenceId,
        "sequenceName": sequenceType + "_" + sequenceName,
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

    sequences[sequenceType + "_" + sequenceName] = sequence


def processSequences():
    sequenceId = 1024

    print("Processing sequences from FASTA")

    with open(FASTA, "r") as fasta:
        while line := fasta.readline().rstrip():
            if line.startswith(">"):
                sequenceId += 1
                processSequence(line, sequenceId)

    print(str(sequenceId - 1024) + " sequences processed from FASTA")


processSequences()
# print(json.dumps(sequences, indent=2,sort_keys=False))
processFeatures()
