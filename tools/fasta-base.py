#!/usr/bin/env python

import sys


def usage():
    print("usage: " + sys.argv[0] + " <fasta> <position>")
    print("    or " + sys.argv[0] + " <fasta> <start> <end>")
    exit(1)


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        usage()

    file_path = sys.argv[1]
    start_position = int(sys.argv[2]) - 1
    end_position = int(sys.argv[3]) if len(sys.argv) > 3 else start_position + 1

    with open(file_path, "r") as fasta:
        contents = "".join(
            [line.strip() for line in fasta.readlines() if not line.strip().startswith(">")]
        )
        print(contents[start_position:end_position])


if __name__ == "__main__":
    main()
