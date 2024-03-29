#!/usr/bin/env python3

from pprint import pprint

from argparse import ArgumentParser, Namespace
import threading
import copy
import gzip
import sys
import os
import json

from pe_read_maker import pe_make_read
from api import find_contig
import genome_operations


def defineArguments():
    parser = ArgumentParser()

    parser.add_argument(
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="Sample name",
    )

    parser.add_argument(
        "--reference",
        required=True,
        action="store",
        metavar="REFERENCE",
        dest="reference",
        help="Full path to reference FNA",
    )

    parser.add_argument(
        "--sequences",
        nargs="+",
        action="store",
        metavar="SEQUENCES",
        dest="sequences",
        default=[],
    )

    parser.add_argument(
        "--read-length",
        action="store",
        metavar="READ_LENGTH",
        dest="readLength",
        default=151,
        type=int,
    )

    parser.add_argument(
        "--read-count",
        action="store",
        metavar="READ_COUNT",
        dest="readCount",
        default=250_000,
    )

    parser.add_argument(
        "--print-queue",
        action="store_true",
        dest="printQueue",
        default=False,
        help="Print operations queue",
    )

    parser.add_argument(
        "--print-reads",
        action="store_false",
        dest="printReads",
        default=True,
        help="Generate R1/R2 files",
    )

    parser.add_argument(
        "--write-fasta",
        action="store_true",
        dest="writeFasta",
        default=False,
        help="Write FASTA files for target sequences",
    )

    parser.add_argument(
        "--mutations",
        nargs="+",
        required=False,
        action="store",
        metavar="MUTATIONS",
        dest="mutations",
        help="JSON files containing list of mutations to apply per sequence",
        default=[],
    )

    return vars(parser.parse_args())


def read_fasta(fa):
    """
    Reads the specified FASTA / multi-FASTA and returns an array of
    contigs with their name, genome, and length.
    """
    result = {}

    with open(fa, "r") as f:
        name = ""

        while line := f.readline().rstrip():
            if line.startswith(">"):
                name = line.removeprefix(">")
                if not name in result:
                    result[name] = []
            else:
                result[name].append(line)

    fasta = [{"contig": k, "accession": k.split(" ")[0], "genome": "".join(v)} for k, v in result.items()]

    return [
        {
            "contig": fa["contig"],
            "accession": fa["accession"],
            "genome": fa["genome"],
            "length": len(fa["genome"]),
        }
        for fa in fasta
    ]


def make_genome(fa, viromes, read_length):
    genome = []

    fa_clone = copy.deepcopy(fa)

    for clone in fa_clone:
        clone["genome"] = list(clone["genome"] + "".join(["N" for _ in range(read_length + 1)]))

        # the list of operations is based on the length of the genome, each
        # position can have at most one variant applied. the length of
        # operations is not affected by the length of the genome during processing
        # because we work backwards
        clone["operations"] = [None for _ in range(clone["length"])]

        genome.append(clone)

    return [g for g in genome if g["accession"] in viromes]


def print_operation(contig, operations):
    print("Reporting " + str(len(operations)) + " for contig " + contig["accession"])
    for op in operations:
        print(op)


def print_operation_queues(genome):
    for index in range(len(genome)):
        contig = genome[index]
        operations = list(reversed([op for ndx, op in enumerate(contig["operations"]) if op != None]))
        if len(operations) > 0:
            print_operation(contig, operations)


def apply_operation_queue(contig, operations, target_read_length):
    print("Applying " + str(len(operations)) + " operations to contig " + contig["accession"])
    for ops in operations:
        for op in ops:
            match (op["verb"]):
                case "ts":
                    genome_operations.do_transition_variant(contig, op["position"])
                case "tv":
                    genome_operations.do_transversion_variant(contig, op["position"], op["slot"])
                case "snv":
                    genome_operations.do_explicit_snv(contig, op["position"], op["base"])
                case "mnv":
                    genome_operations.do_explicit_mnv(contig, op["position"], op["bases"])
                case "ins":
                    genome_operations.do_insert(contig, op["position"], op["bases"])
                case "del":
                    genome_operations.do_delete(contig, op["position"], op["length"])

    # now we correct the length because of ins / del operations
    # and we need to remove the N tail
    contig["length"] = len(contig["genome"]) - target_read_length


def apply_operation_queues(genome, target_read_length):
    for index in range(len(genome)):
        contig = genome[index]
        operations = list(reversed([op for ndx, op in enumerate(contig["operations"]) if op != None]))
        if len(operations) > 0:
            apply_operation_queue(contig, operations, target_read_length)


def apply_mutations(contig, mutations, target_read_length):
    operations = list(reversed([op for ndx, op in enumerate(mutations) if op != None]))
    if len(operations) > 0:
        apply_operation_queue(contig, operations, target_read_length)


def print_reads(sample_name, reads, read_count, direction):
    print("Writing R" + str(direction))
    f = gzip.open(sample_name + "_R" + str(direction) + ".fastq.gz", "wb", compresslevel=6)
    for r in range(read_count):
        read = reads[r]["r" + str(direction)]

        f.write(
            bytes(
                read["information"] + "_" + str(r + 1) + " " + str(direction) + ":N:0:1\n",
                "ascii",
            )  # type: ignore
        )
        f.write(bytes(read["read"] + "\n", "ascii"))  # type: ignore
        f.write(bytes("+\n", "ascii"))  # type: ignore
        f.write(bytes(read["quality"] + "\n", "ascii"))  # type: ignore
    f.close()


def write_fasta(filename, sample_name, fasta):
    with open(filename, "w") as f:
        if sample_name.endswith("!"):
            f.write(">" + str(sample_name.removesuffix("!")) + "\n")
        else:
            f.write(">" + str(sample_name) + " | (" + fasta["accession"] + ")\n")

        f.writelines([l + "\n" for l in fasta])


def write_genome_files(genome, viromes):
    if len(viromes) == 0:
        print("No viromes selected to output")
        quit(1)

    for index in range(len(genome)):
        contig = genome[index]
        if contig["accession"] in viromes:
            filename = os.path.join(root_folder, contig["accession"]) + ".fa"
            fasta = contig["genome"]

            # strip trailing N (we might have put them there)
            for b in range(len(fasta) - 1, 0, -1):
                if fasta[b] == "N":
                    del fasta[b]

            lines = ["".join(fasta[i : i + 60]) for i in range(0, len(fasta), 60)]
            write_fasta(filename, contig["accession"] + "!", lines)


def make_and_print_reads(read_count, sample, fa, viromes):
    if len(viromes) == 0:
        print("No viromes selected to construct")
        quit(1)

    print("Allocating read array for " + str(sample))
    reads = [{}] * read_count

    print("Building reads")
    for r in range(read_count):
        reads[r] = pe_make_read(fa, target_read_length, viromes)

    threads = []

    for r in range(2):
        x = threading.Thread(
            target=print_reads,
            args=(
                str(sample),
                reads,
                read_count,
                r + 1,
            ),
        )
        threads.append(x)
        x.start()

    for _, thread in enumerate(threads):
        thread.join()


def dump_operations(contig, sequence):
    mutations = {sequence: [x for x in contig["operations"] if x != None]}
    print(json.dumps(mutations, indent=2))


if __name__ == "__main__":
    root_folder = os.path.dirname(sys.argv[0])

    options = defineArguments()

    target_read_length = int(options["readLength"])

    # load and create the sample's genome
    genome = make_genome(
        read_fasta(options["reference"]),
        options["sequences"],
        target_read_length,
    )

    if options["mutations"] != None:
        for mutation_file in options["mutations"]:
            with open(mutation_file) as f:
                mutations = json.load(f)
                for sequence in mutations.keys():
                    contig = find_contig(genome, sequence)
                    apply_mutations(contig, mutations[sequence], target_read_length)

    if options["printQueue"]:
        print_operation_queues(genome)

    apply_operation_queues(genome, target_read_length)

    if options["printReads"]:
        make_and_print_reads(int(options["readCount"]), options["sample"], genome, options["sequences"])

    if options["writeFasta"]:
        write_genome_files(genome, options["sequences"])
