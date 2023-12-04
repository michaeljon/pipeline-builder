#!/usr/bin/env python

from colorama import Fore
import csv
import gzip
import re
import sys


def usage():
    print("usage: " + sys.argv[0] + " <vcf> <vcf>")
    exit(1)


def read(vcf):
    fieldnames = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    result = {}

    if vcf.endswith("gz"):
        with gzip.open(vcf, "rt") as f:
            r = csv.DictReader(
                filter(lambda row: str(row[0]) != "#", f),  # type: ignore
                fieldnames,
                delimiter="\t",
            )
            for variants in r:
                result[int(variants["POS"])] = variants
    else:
        with open(vcf, "rt") as f:
            r = csv.DictReader(
                filter(lambda row: str(row[0]) != "#", f),  # type: ignore
                fieldnames,
                delimiter="\t",
            )
            for variants in r:
                result[int(variants["POS"])] = {
                    "CHROM": variants["CHROM"],
                    "POS": variants["POS"],
                    "REF": variants["REF"],
                    "ALT": variants["ALT"],
                }

    return result


def is_match(left, right):
    return left["REF"] == right["REF"] and left["ALT"] == right["ALT"]


def get_difference(left, right):
    delta = "{0:<10}\t{1:<6}\t{2:<5}\t{3:<5}\t{4:<10}\t{5:<6}\t{6:<5}\t{7:<5}".format(
        left["CHROM"],
        left["POS"],
        left["REF"],
        left["ALT"],
        right["CHROM"],
        right["POS"],
        right["REF"],
        right["ALT"],
    )

    return delta

    if left["POS"] != "" and right["POS"] != "":
        return Fore.RED + delta + Fore.RESET
    else:
        return delta


def get_contig_from_line(line):
    pattern = re.compile("##contig=<ID=([^,]+),length=([0-9]+)>")
    m = pattern.match(line)
    if m != None:
        return m.group(1)


def get_contig_from_vcf(vcf):
    if vcf.endswith("gz"):
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if str(line).startswith("##contig"):
                    return get_contig_from_line(line)
    else:
        with open(vcf, "rt") as f:
            for line in f:
                if line.startswith("##contig"):
                    return get_contig_from_line(line)


def compare(left, left_contig, right, right_contig):
    differences = {}

    dummy = {
        "CHROM": "",
        "POS": "",
        "REF": "",
        "ALT": "",
    }

    left_dummy = {**dummy, "CHROM": left_contig}
    right_dummy = {**dummy, "CHROM": right_contig}

    # look for matching positions where the ref/alt pair is different
    for pos in left.keys():
        if pos in right:
            # both have an entry
            if not is_match(left[pos], right[pos]):
                differences[int(pos)] = get_difference(left[pos], right[pos])

    # look for variants in the left vcf that aren't in the right
    for pos in left.keys():
        if not pos in right:
            # left-only variant
            differences[int(pos)] = get_difference(left[pos], right_dummy)

    # look for variants in the right vcf that aren't in the left
    for pos in right.keys():
        if not pos in left:
            # right-only variant
            differences[int(pos)] = get_difference(left_dummy, right[pos])

    print(
        "{0:<10}\t{1:<6}\t{2:<5}\t{3:<5}\t{4:<10}\t{5:<6}\t{6:<5}\t{7:<5}".format(
            "CHROM   ", "POS", "REF", "ALT", "CHROM   ", "POS", "REF", "ALT"
        )
    )

    for difference in dict(sorted(differences.items())):
        print(differences[difference])


def main():
    if len(sys.argv) != 3:
        usage()

    left = read(sys.argv[1])
    left_contig = get_contig_from_vcf(sys.argv[1])
    right = read(sys.argv[2])
    right_contig = get_contig_from_vcf(sys.argv[2])

    compare(left, left_contig, right, right_contig)


if __name__ == "__main__":
    main()
