#!/usr/bin/env python

import sys
import json
import gzip
import numpy as np

if len(sys.argv) != 4:
    print("usage: " + sys.argv[0] + " <sample> <organism> <pileup>")
    exit(1)

feature_data = {}
with open("hcov-regions.json") as f:
    feature_data = json.load(f)

sample = sys.argv[1]
organism = sys.argv[2]
pileup = sys.argv[3]

if organism not in feature_data:
    print("Unknown organism " + organism + " check region data")
    exit(2)

sequence_id = feature_data[organism]["sequence"]
genes = feature_data[organism]["genes"]
region_data = {k: [0 for _ in range(v["start"], v["stop"] + 1)] for i, (k, v) in enumerate(genes.items())}

with gzip.open(pileup, "rt") as f:
    while line := f.readline().rstrip():
        if line.startswith("#") == False:
            elements = line.split("\t")

            if elements[0] == sequence_id:
                pos = int(elements[1])
                infos = elements[7].split(";")
                for i in infos:
                    v = i.split("=")
                    if i.startswith("DP="):
                        for g in genes:
                            if genes[g]["start"] <= pos and pos <= genes[g]["stop"]:
                                region_data[g][pos - genes[g]["start"]] = int(v[1])

print("sample\torganism\tgene\tmean\tmedian\tmin\tmax\tstdev")
for gene in genes:
    mean = np.mean(region_data[gene])
    median = np.median(region_data[gene])
    min = np.min(region_data[gene])
    max = np.max(region_data[gene])
    stdev = np.std(region_data[gene])

    print(
        "{}\t{}\t{}\t{:.1f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.4f}".format(
            sample, organism, gene, mean, median, min, max, stdev
        )
    )
