from time import time
from random import random, randint

from rules import complements


def pe_make_read(fa, read_len):
    contig = randint(0, len(fa) - 1)

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
    return "".join(["I" for _ in range(len(read))])
