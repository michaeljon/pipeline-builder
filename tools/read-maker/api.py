# an operation is
#   location: contig, position
#   verb:     replace-one + base
#             replace-many + [base]
#             ins + [bases]
#             del + count
#             indel / delins + count + [bases]

# an operation list is
#  order list (descending) by position of operation


def find_contig(genome, accession):
    for index in range(len(genome)):
        contig = genome[index]
        if contig["accession"] == accession:
            return contig

    print("Attempt to cause variant in accession " + accession)
    quit(1)


def queue_transition_variant(genome, accession, position):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append({"verb": "ts", "position": position})


def queue_tranversion_variant(genome, accession, position, slot):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append({"verb": "tv", "position": position, "slot": slot})


def queue_explicit_snv(genome, accession, position, base):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append({"verb": "snv", "position": position, "base": base})


def queue_explicit_mnv(genome, accession, position, bases):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append(
        {
            "verb": "mnv",
            "position": position,
            "bases": bases,
        }
    )


def queue_insert(genome, accession, position, bases):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append(
        {
            "verb": "ins",
            "position": position,
            "bases": bases,
        }
    )


def queue_delete(genome, accession, position, length):
    contig = find_contig(genome, accession)
    if contig["operations"][position] == None:
        contig["operations"][position] = []

    contig["operations"][position].append(
        {
            "verb": "del",
            "position": position,
            "length": length,
        }
    )
