#!/usr/bin/env python

import sys


def usage():
    print("usage: " + sys.argv[0] + " <fasta> <sequence>")
    exit(1)


def main():
    if len(sys.argv) != 3:
        usage()

    file_path = sys.argv[1]
    sequence = sys.argv[2]

    with open(file_path, "r") as fasta:
        contents = "".join([line.strip() for line in fasta.readlines() if not line.strip().startswith(">")])

        locus = contents.find(sequence)
        while locus >= 0:
            print(locus)
            locus = contents.find(sequence, locus + 1)


if __name__ == "__main__":
    main()
