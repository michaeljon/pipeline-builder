#!/usr/bin/env python3

from pprint import pprint

import threading
import copy
import gzip

from pe_read_maker import pe_make_read
import genome_operations

target_read_length = 151


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
        {"contig": fa["contig"], "accession": fa["accession"], "genome": fa["genome"], "length": len(fa["genome"])}
        for fa in fasta
    ]


def make_genome(fa, read_length):
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

    return genome


def find_contig(genome, accession):
    for index in range(len(genome)):
        contig = genome[index]
        if contig["accession"] == accession:
            return contig

    print("Attempt to cause variant in accession " + accession)
    quit(1)


# an operation is
#   location: contig, position
#   verb:     replace-one + base
#             replace-many + [base]
#             ins + [bases]
#             del + count
#             indel / delins + count + [bases]

# an operation list is
#  order list (descending) by position of operation


def queue_transition_variant(genome, accession, position):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "ts", "position": position}


def queue_tranversion_variant(genome, accession, position, slot):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "tv", "position": position, "slot": slot}


def queue_explicit_snv(genome, accession, position, base):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "snv", "position": position, "base": base}


def queue_explicit_mnv(genome, accession, position, bases):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "mnv", "position": position, "bases": bases}


def queue_insert(genome, accession, position, bases):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "ins", "position": position, "bases": bases}


def queue_delete(genome, accession, position, length):
    contig = find_contig(genome, accession)
    contig["operations"][position] = {"verb": "del", "position": position, "length": length}


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


def apply_operation_queue(contig, operations):
    print("Applying " + str(len(operations)) + " operations to contig " + contig["accession"])
    for op in operations:
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


def apply_operation_queues(genome):
    for index in range(len(genome)):
        contig = genome[index]
        operations = list(reversed([op for ndx, op in enumerate(contig["operations"]) if op != None]))
        if len(operations) > 0:
            apply_operation_queue(contig, operations)


def print_reads(reads, read_count, direction):
    print("Writing R" + str(direction))
    f = gzip.open("reads_R" + str(direction) + ".fastq.gz", "wb", compresslevel=6)
    for r in range(read_count):
        read = reads[r]["r" + str(direction)]

        f.write(bytes(read["information"] + "_" + str(r + 1) + "/" + str(direction) + "\n", "ascii"))
        f.write(bytes(read["read"] + "\n", "ascii"))
        f.write(bytes("+\n", "ascii"))
        f.write(bytes(read["quality"] + "\n", "ascii"))
    f.close()


def print_reads_debug(reads, read_count, direction):
    print("Writing R" + str(direction))
    f = open("/dev/stdout", "w")
    for r in range(read_count):
        read = reads[r]["r" + str(direction)]

        f.write(read["information"] + "_" + str(r + 1) + "/" + str(direction) + "\n")
        f.write(read["read"] + "\n")
        f.write("+\n")
        f.write(read["quality"] + "\n")
    f.close()


def make_and_print_reads(fa):
    read_count = 1_000_000

    print("Allocating read array")
    reads = [{}] * read_count

    print("Building reads")
    for r in range(read_count):
        reads[r] = pe_make_read(fa, target_read_length)

    threads = []

    for r in range(2):
        x = threading.Thread(
            target=print_reads,
            args=(
                reads,
                read_count,
                r + 1,
            ),
        )
        threads.append(x)
        x.start()

    for _, thread in enumerate(threads):
        thread.join()


# load and create the sample's genome
genome = make_genome(
    read_fasta("/Users/michaeljon/reference/hcov-panel/hcov-panel.fna"),
    target_read_length,
)

# # nl63
# contig = find_contig(genome, "AY567487.2")
# for n in range(1000, contig["length"], 1500):
#     queue_transition_variant(genome, "AY567487.2", n)

# # oc43
# queue_delete(genome, "AY585228.1", 100, 10)
# queue_explicit_mnv(genome, "AY585228.1", 4001, list("ATCGATCGATCGATCGATCG"))
# queue_transition_variant(genome, "AY585228.1", 3995)
# queue_transition_variant(genome, "AY585228.1", 3997)
# queue_transition_variant(genome, "AY585228.1", 3999)
# queue_delete(genome, "AY585228.1", 100, 10)

# 229e
contig = find_contig(genome, "AF304460.1")
for n in range(24500, 24750, 25):
    queue_transition_variant(genome, "AF304460.1", n)
    queue_delete(genome, "AF304460.1", n + 10, 2)
queue_insert(genome, "AF304460.1", 24200, list("AAAATTTTCCCCGGGGATCGATCGATCGATCG"))

# # hku1
# queue_tranversion_variant(genome, "AY597011.2", 3995, 0)
# queue_tranversion_variant(genome, "AY597011.2", 3997, 0)
# queue_tranversion_variant(genome, "AY597011.2", 3999, 0)

# print_operation_queues(genome)
apply_operation_queues(genome)
make_and_print_reads(genome)
