from time import time
from random import random, randint

from rules import complements


def pe_make_read(fa, read_len, viromes):
    contig = randint(0, len(viromes) - 1)

    forward = perfect_read(fa[contig], read_len)
    reverse = complement(reversed(forward))
    information = "@" + fa[contig]["accession"] + "_" + str(int(time()))

    return {
        "contig": fa[contig]["accession"],
        "r1": {"information": information, "read": forward, "quality": quality(forward)},
        "r2": {"information": information, "read": reverse, "quality": quality(forward)},
    }


def perfect_read(contig, read_len):
    position = randint(0, contig["length"] - 1)
    read = "".join(contig["genome"][position : position + read_len])

    if len(read) != read_len:
        print("Attempted to read beyond end of genome")
        quit(1)

    return read


def complement(read):
    return "".join([complements[base] for base in list(read)])


def quality(read):
    # this is a phred33 score of 30, it's the highest non-overlapping and
    # unambigious score possible, see
    # https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html
    return "".join(["?" for _ in range(len(read))])
