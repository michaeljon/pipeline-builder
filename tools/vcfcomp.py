from colorama import Fore
import csv
import gzip
import sys


def usage():
    print("usage: " + sys.argv[0] + " <vcf> <vcf>")
    exit(1)

def read(vcf):
    fieldnames=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]    
    result = {}

    if vcf.endswith("gz"):
        with gzip.open(vcf, "rt") as f:
            r = csv.DictReader(filter(lambda row: str(row[0]) != "#", f),  # type: ignore
                            fieldnames,
                            delimiter="\t",)
            for variants in r:
                result[int(variants["POS"])] = variants
    else:
        with open(vcf, "rt") as f:
            r = csv.DictReader(filter(lambda row: str(row[0]) != "#", f),  # type: ignore
                            fieldnames,
                            delimiter="\t",)
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
    delta = "{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}".format(
          left["CHROM"], left["POS"], left["REF"], left["ALT"],
          right["CHROM"], right["POS"], right["REF"], right["ALT"],)
    
    if left["POS"] != "" and right["POS"] != "":
        return Fore.RED + delta + Fore.RESET
    else:
        return delta

def compare(left, right):
    dummy = {
                    "CHROM": "",
                    "POS": "",
                    "REF": "",
                    "ALT": "",
                    }
    differences = {}

    left_chrom = next(iter(left.values()))["CHROM"]
    right_chrom = next(iter(right.values()))["CHROM"]


    for pos in left.keys():
        if pos in right:
            # both have an entry
            if not is_match(left[pos], right[pos]):
                differences[int(pos)] = get_difference(left[pos], right[pos])

    for pos in left.keys():
        if not pos in right:
            # left-only variant
            differences[int(pos)] = get_difference(left[pos], {**dummy, "CHROM": right_chrom})
            

    for pos in right.keys():
        if not pos in left:
            # right-only variant
            differences[int(pos)] = get_difference({**dummy, "CHROM": left_chrom}, right[pos])

    print("{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}".format(
          "CHROM   ", "POS", "REF", "ALT",
          "CHROM   ", "POS", "REF", "ALT"))
    
    for difference in dict(sorted(differences.items())):
        print(differences[difference])


def main():
    if len(sys.argv) != 3:
        usage()

    left = read(sys.argv[1])
    right = read(sys.argv[2])

    compare(left, right)



if __name__ == "__main__":
    main()
