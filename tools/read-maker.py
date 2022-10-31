#!/usr/bin/env python

from random import random, randint
import gzip

from matplotlib.font_manager import json_dump

reference = []
genome = []

adapters = [
    "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",  # Illumina Paired End Adapter 1
    "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",  # Illumina Paired End Adapter 2
    # "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",  # Illumina Paired End Adapter 1
    # "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",  # Illumina Paired End Adapter 2
    # "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",  # Illumina Paried End PCR Primer 1
    # "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",  # Illumina Paired End PCR Primer 2
    # "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",  # Illumina Paried End Sequencing Primer 1
    # "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",  # Illumina Paired End Sequencing Primer 2
]

transpositions = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
atcgn = ["A", "T", "C", "G", "N"]
atcg = ["A", "T", "C", "G"]

unmapped = 0
sequencer_variations = 0
transpose_variations = 0


def readReference(ref):
    with open(ref, "r") as f:
        for l in f.readlines():
            if not l.startswith(">"):
                reference.extend(l.strip())
                genome.extend(l.strip())


def createSampleGenome(dogc: bool, dongc: bool, domnv: bool, doins: bool, dodels: bool):
    global genome

    gc = 0
    nongc = 0
    mnv = 0
    insertions = 0
    deletions = 0

    indels = []

    # step 1 - gc-conserve randomize 1/1000 of the sequence
    if dogc == True:
        for n in range(0, len(genome) - 1):
            if random() < 0.001:
                gc += 1
                indels.append(
                    {"type": "gc", "length": 1, "position": n + 1, "from": genome[n], "to": transpositions[genome[n]]}
                )
                genome[n] = transpositions[genome[n]]

    # step 2 - randomize 1/1000 of the sequence
    if dongc == True:
        for n in range(0, len(genome) - 1):
            if random() < 0.001:
                nongc += 1
                atgc = [ch for ch in atcg if ch != genome[n]]
                ch = atgc[randint(0, len(atgc) - 1)]
                indels.append({"type": "ngc", "length": 1, "position": n + 1, "from": genome[n], "to": ch})
                genome[n] = ch

    # step 3 - perform an MNV approx 3/10000
    if domnv == True:
        for n in range(0, len(genome) - 1):
            if random() < 0.0003:
                mnv += 1
                l = randint(2, 5)
                ins = [atcg[randint(0, len(atcg) - 1)] for _ in range(0, l)]

                indels.append({"type": "mnv", "length": l, "position": n + 1, "from": genome[n : n + l], "to": ins})

                for p in range(min(len(ins), len(genome)) - 1):
                    genome[n + p] = ins[p]

    # step 4 - insert random number of sequences of random lengths
    if doins == True:
        toapply = randint(5, 10)
        for n in range(toapply):
            insertions += 1
            l = randint(3, 10)
            ins = [atcg[randint(0, len(atcg) - 1)] for _ in range(0, l)]
            pos = randint(0, len(genome) - 1)

            indels.append({"type": "insertion", "length": l, "ins": ins, "position": pos + 1})

            tmp = genome[0:pos]
            tmp += ins
            tmp += genome[pos:]

            genome = tmp

    # step 5 - delete random number of sequences of random lengths
    if dodels == True:
        toapply = randint(5, 10)
        for n in range(toapply):
            deletions += 1
            l = randint(3, 10)
            pos = randint(0, len(genome) - 1)

            indels.append({"type": "deletion", "length": l, "position": pos + 1})

            tmp = genome[0:pos]
            tmp += genome[pos + l + 1 :]

            genome = tmp

    if len(genome) < len(reference):
        genome.extend(["A" for _ in range(len(reference) - len(genome))])

    print("Stats")
    print("  Reference len " + str(len(reference)))
    print("  Genome len " + str(len(genome)))
    print("  GC-conserving " + str(gc))
    for v in sorted(indels, key=lambda i: i["position"]):
        if v["type"] == "gc":
            print("   " + str(v["position"]) + ": " + v["from"] + " -> " + v["to"])

    print("  Non-conserving " + str(nongc))
    for v in sorted(indels, key=lambda i: i["position"]):
        if v["type"] == "ngc":
            print("   " + str(v["position"]) + ": " + v["from"] + " -> " + v["to"])

    print("  MNV " + str(mnv))
    for v in sorted(indels, key=lambda i: i["position"]):
        if v["type"] == "mnv":
            print("   " + str(v["position"]) + ": " + "".join(v["from"]) + " -> " + "".join(v["to"]))

    print("  Insertions " + str(insertions))
    for v in sorted(indels, key=lambda i: i["position"]):
        if v["type"] == "insertion":
            print("   " + str(v["position"]) + ": " + "".join(v["ins"]))

    print("  Deletions " + str(deletions))
    for v in sorted(indels, key=lambda i: i["position"]):
        if v["type"] == "deletion":
            print("   " + str(v["position"]) + ": " + str(v["length"]))


def writeFasta(fasta):
    fasta.write(">SAMPLE (ref SYN)\n")

    for chunk in [genome[i : i + 80] for i in range(0, len(genome), 80)]:
        fasta.write("".join(chunk) + "\n")


def makeName(instrument, run, flowcell, lane, tile, x, y, member) -> str:
    # @A01343:127:H5CHVDSX5:3:1101:3260:1157 1:N:0:CAACAATG+CCGTGAAG
    return (
        "@"
        + instrument
        + ":"
        + str(run)
        + ":"
        + flowcell
        + ":"
        + str(lane)
        + ":"
        + str(tile)
        + ":"
        + str(x)
        + ":"
        + str(y)
        + " "
        + str(member)
        + ":N:0:CAACAATG+CCGTGAAG"
    )


def selectFragment(length):
    global sequencer_variations
    global unmapped

    # about 75/1000 generate a completely unmapped read (this means about
    # 5% of our genetic source material isn't part of our organism)
    if random() < 0.075:
        unmapped += 1
        return "".join([atcg[randint(0, len(atcg) - 1)] for _ in range(length)])

    start = randint(0, len(genome))
    read = genome[start : start + length]

    # if we read at the tail end of the genome we'll return an unmapped read instead
    if len(read) == 0:
        unmapped += 1
        return "".join([atcg[randint(0, len(atcg) - 1)] for _ in range(length)])

    # about 5% of the time we will have an "issue" with the current read
    if random() < 0.05:
        for n in range(0, len(read) - 1):
            # and we'll walk over the entire read and about 2% of the
            # time we'll randomly change a read
            if random() < 0.02:
                sequencer_variations += 1
                read[n] = atcgn[randint(0, len(atcgn) - 1)]

    # about 1% of the time we'll attach a polyX read between 10 and 20 bases long
    if random() < 0.01:
        polyX = "G"
        read += [polyX for _ in range(randint(10, 20))]

    return "".join(read)


def clamp(n, smallest, largest):
    return max(smallest, min(n, largest))


def makeQuality(read):
    length = len(read)
    quality = []

    for ch in range(length):
        if read[ch] == "N":
            quality.append("!")
        else:
            quality.append(chr(randint(ord(":"), ord("I"))))

    return "".join(quality)


def transpose(read):
    global transpose_variations

    transposed_read = []

    for ch in read:
        # assume it's clean
        n = transpositions[ch]

        # wipe it out otherwise
        if random() < 0.0001:
            transpose_variations += 1
            ch = atcgn[randint(0, len(atcgn) - 1)]

        transposed_read.append(n)

    return "".join(transposed_read)


def writeRead(f, name, read, quality):
    f.write(bytes(name + "\n" + read + "\n+\n" + quality + "\n", "ascii"))


def makeTile():
    return 1011


def makeX():
    return randint(1111, 9999)


def makeY():
    return randint(1111, 9999)


def makePixelDistance():
    return randint(-10, 10)


r1 = gzip.open("/Users/michaeljon/tmp/fastq/SAMPLE_R1_001.fastq.gz", "wb", compresslevel=6)
r2 = gzip.open("/Users/michaeljon/tmp/fastq/SAMPLE_R2_001.fastq.gz", "wb", compresslevel=6)

# r1 = open("/Users/michaeljon/tmp/fastq/SAMPLE_R1_001.fastq", "wb")
# r2 = open("/Users/michaeljon/tmp/fastq/SAMPLE_R2_001.fastq", "wb")

readReference("/Volumes/Genomics/pipeline/reference/synthetic/synthetic.fna")
createSampleGenome(dogc=True, dongc=True, domnv=True, doins=True, dodels=True)
with open("/Users/michaeljon/tmp/fastq/SAMPLE.fasta", "w") as fasta:
    writeFasta(fasta)

inst = "MN87AB"
run = randint(1001, 2001)
flowcell = "H5CHVDSX5"
lane = randint(1001, 2001)

bases = 0
gcConserving = 0
randomSnv = 0
insertSv = 0
deletionSv = 0
adapters_added = 0

reads = 0
non_optical_duplicates = 0
optical_duplicates = 0

for n in range(50_000):  # range(randint(50_000, 100_000)):
    readlength = 151
    add_adapter = random() < 0.1
    if add_adapter == True:
        adapters_added += 1
        readlength -= len(adapters[0])

    reads += 1
    bases += readlength

    tile = makeTile()
    x = makeX()
    y = makeX()

    fragment = selectFragment(readlength)
    read = fragment
    if add_adapter == True:
        read = fragment + adapters[0]

    forward_name = makeName(inst, run, flowcell, lane, tile, x, y, 1)
    forward_read = read
    forward_quality = makeQuality(forward_read)

    read = fragment
    if add_adapter == True:
        read = transpose(fragment[::-1]) + adapters[1]

    reverse_name = makeName(inst, run, flowcell, lane, tile, x, y, 2)
    reverse_read = read
    reverse_quality = makeQuality(reverse_read)

    writeRead(r1, forward_name, forward_read, forward_quality)
    writeRead(r2, reverse_name, reverse_read, reverse_quality)

    # duplicates (optical)
    # about 1/100 are sequencer issues
    if random() < 0.01:
        bases += readlength

        optical_duplicates += 1

        dist = makePixelDistance()

        forward_name = makeName(inst, run, flowcell, lane, tile, x + dist, y + dist, 1)
        reverse_name = makeName(inst, run, flowcell, lane, tile, x + dist, y + dist, 2)

        writeRead(r1, forward_name, forward_read, forward_quality)
        writeRead(r2, reverse_name, reverse_read, reverse_quality)

r2.close()
r1.close()

print("Reads " + str(reads))
print("Bases " + str(bases))
print("Adapters added " + str(adapters_added))
# print("Non-optical duplicates " + str(non_optical_duplicates))
print("Optical duplicates " + str(optical_duplicates))
print("Unmapped reads " + str(unmapped))
print("Sequencer variations applied " + str(sequencer_variations))
print("Transpose variations applied " + str(transpose_variations))
# print("L/W depth " + str(int((reads * readlength * 2) / len(reference))))
