#!/usr/bin/env python

import sys
import os


def read_fasta(fa):
    split_size = 75

    with open(fa, "r") as f:
        lines = []
        name = ""

        while line := f.readline().rstrip():
            if line.startswith(">"):
                name = line.removeprefix(">")
            else:
                lines.append(line)

        fasta = "".join(lines)
        lines = [fasta[i : i + split_size] for i in range(0, len(fasta), split_size)]

        return {"name": name, "fasta": lines}


def write_fasta(filename, sample_name, fasta):
    with open(filename, "w") as f:
        if sample_name.endswith("!"):
            f.write(">" + str(sample_name.removesuffix("!")) + "\n")
        else:
            f.write(">" + str(sample_name) + " | (" + fasta["name"] + ")\n")

        f.writelines([l + "\n" for l in fasta["fasta"]])


def usage():
    print("usage: " + sys.argv[0] + " <fa> [label]")
    print("example:\n\t" + sys.argv[0] + "sample-01-hku1.fa sample-01-hku1")
    exit(1)


filename = None
sample_name = None

if len(sys.argv) == 3:
    filename = sys.argv[1]
    sample_name = sys.argv[2]

elif len(sys.argv) == 2:
    filename = sys.argv[1]
    sample_name = os.path.basename(filename.removesuffix(".consensus.fa"))

else:
    usage()

write_fasta(filename, sample_name, read_fasta(filename))
