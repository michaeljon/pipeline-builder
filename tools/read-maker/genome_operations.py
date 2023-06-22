from rules import transitions, transversions


def do_transition_variant(contig, position):
    current_base = contig["genome"][position - 1]
    contig["genome"][position - 1] = transitions[current_base]

    print(
        "{}: ts @{} from {} to {}".format(
            contig["accession"],
            position - 1,
            current_base,
            contig["genome"][position - 1],
        )
    )


def do_transversion_variant(contig, position, slot):
    current_base = contig["genome"][position - 1]
    contig["genome"][position - 1] = transversions[current_base][slot]

    print(
        "{}: tv @{} from {} to {}".format(
            contig["accession"],
            position - 1,
            current_base,
            contig["genome"][position - 1],
        )
    )


def do_explicit_snv(contig, position, base):
    current_base = contig["genome"][position - 1]
    contig["genome"][position - 1] = base

    print(
        "{}: snv @{} from {} to {}".format(
            contig["accession"],
            position - 1,
            current_base,
            transitions[current_base],
        )
    )


def do_explicit_mnv(contig, position, bases):
    current_bases = contig["genome"][position - 1 : position - 1 + len(bases) - 1]
    contig["genome"][position - 1 : position - 1 + len(bases) - 1] = bases

    print(
        "{}: mnv @{} from {} to {}".format(
            contig["accession"],
            position - 1,
            "[" + ",".join(current_bases) + "]",
            "[" + ",".join(bases) + "]",
        )
    )


def do_insert(contig, position, bases):
    print(
        "{}: ins @{} of {}".format(
            contig["accession"],
            position,
            "[" + ",".join(bases) + "]",
        )
    )

    tmp = contig["genome"][0 : position - 1]
    tmp += bases
    tmp += contig["genome"][position - 1 :]

    contig["genome"] = tmp


def do_delete(contig, position, count):
    print(
        "{}: del @{} of {} '{}'".format(
            contig["accession"],
            position,
            count,
            contig["genome"][position - 1 : position - 1 + count],
        )
    )

    del contig["genome"][position - 1 : position - 1 + count]
