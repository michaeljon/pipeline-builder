#!/usr/bin/env python

import argparse
import csv
import gzip

from os.path import exists, expandvars
from sys import stderr

parser = argparse.ArgumentParser()
parser.add_argument(
    "-d",
    "--dbsnp",
    required=True,
    action="store",
    metavar="VCF",
    dest="vcf",
    help="Path to dbsnp VCF",
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
parser.add_argument(
    "-l",
    "--limit",
    required=False,
    action="store",
    metavar="LIMIT",
    dest="limit",
    default=0,
    type=int,
    help="Limit number of input variants to process",
)
parser.add_argument(
    "-c",
    "--rare-variants",
    required=False,
    action="store_true",
    dest="rare-variants",
    default=False,
    help="Process rare variants",
)
opts = parser.parse_args()
options = vars(opts)

if options["tsv"] == "-":
    options["tsv"] = "/dev/stdout"

for opt in ["vcf", "tsv"]:
    options[opt] = expandvars(options[opt])

if exists(options["vcf"]) == False:
    print("Unable to find VCF (--input) at {PATH}".format(PATH=options["vcf"]))
    quit(1)

lookup = [
    "NC_000001.11",
    "NC_000002.12",
    "NC_000003.12",
    "NC_000004.12",
    "NC_000005.10",
    "NC_000006.12",
    "NC_000007.14",
    "NC_000008.11",
    "NC_000009.12",
    "NC_000010.11",
    "NC_000011.10",
    "NC_000012.12",
    "NC_000013.11",
    "NC_000014.9",
    "NC_000015.10",
    "NC_000016.10",
    "NC_000017.11",
    "NC_000018.10",
    "NC_000019.10",
    "NC_000020.11",
    "NC_000021.9",
    "NC_000022.11",
    "NC_000023.11",
    "NC_000024.10",
    "NC_012920.1",
]

##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Variant names from HGVS.    The order of these variants corresponds to the order of the info in the other clinical  INFO tags.">
##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|).  Does not include pseudogenes.">
##INFO=<ID=PSEUDOGENEINFO,Number=1,Type=String,Description="Pairs each of pseudogene symbol:gene id.  The pseudogene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=VC,Number=1,Type=String,Description="Variation Class">
##INFO=<ID=COMMON,Number=0,Type=Flag,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">

# TODO add CLNORIGIN - number / bitmask
# TODO add CLNSIG - number
# TODO add impact tags NSF, NSM, NSN, SYN - boolean (FxnClass)
# TODO add position tags U3, U5, ASS, DSS, INT, R3, R5 - boolean (FxnCode)

print("Processing variants from VCF")

reads_processed = 0

with open(options["tsv"], "w") as tsv:
    writer = csv.writer(tsv, delimiter="\t")
    writer.writerow(
        [
            "reference",
            "sequence_id",
            "position",
            "ref",
            "alt",
            "hgvs_root",
            "hgvs",
            "variant_type",
            "common",
            "gene_info",
            "rs",
        ]
    )

    rows = []

    vcf = gzip.open(options["vcf"], "rt") if options["vcf"].endswith("gz") else open(options["vcf"], "r")

    while line := vcf.readline().rstrip():
        if options["limit"] != 0 and reads_processed > options["limit"]:
            break

        if line.startswith("#") == False:
            reads_processed += 1
            columns = line.split("\t")

            # if not columns[0] in lookup:
            #     continue

            sequence_id = columns[0]
            position = int(columns[1])
            rs = columns[2]
            if len(rs) > 0 and rs.startswith("rs"):
                rs = rs.removeprefix("rs")
            ref = columns[3]
            alt = columns[4]

            info = columns[7].split(";")
            infos = {k: v for k, v in zip([k.split("=")[0] for k in info], [":".join(k.split("=")[1:]) for k in info])}

            common = "1" if "COMMON" in infos else "0"
            if common == "0" and options["rare-variants"] == False:
                continue

            hgvs = infos["CLNHGVS"] if "CLNHGVS" in infos else None
            hgvs_root = hgvs[0 : hgvs.find(",") - 1] if hgvs is not None and "," in hgvs else None

            # yeah, this could be handled better...
            if hgvs is None:
                hgvs = sequence_id + ":g." + str(position)

            if hgvs_root is None:
                hgvs_root = sequence_id + ":g." + str(position)

            variant_type = infos["VC"].lower() if "VC" in infos else None

            # TODO - split this
            gene_info = infos["GENEINFO"] if "GENEINFO" in infos else None
            if gene_info == None:
                gene_info = infos["PSEUDOGENEINFO"] if "PSEUDOGENEINFO" in infos else None

            if rs == "." or rs is None:
                rs = infos["RS"].split("|")[0] if "RS" in infos else None

            rows.append(
                [
                    options["ref"],
                    sequence_id,
                    position,
                    ref,
                    alt,
                    hgvs_root,
                    hgvs,
                    variant_type,
                    common,
                    gene_info,
                    rs,
                ]
            )

            if reads_processed % 10000 == 0:
                writer.writerows(rows)
                rows = []
                stderr.write(".")
                stderr.flush()

    writer.writerows(rows)
    rows = []

    vcf.close()

stderr.write("\n")
print("VCF processed " + str(reads_processed) + " reads")
