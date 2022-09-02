#!/usr/bin/env python

from importlib.resources import read_binary
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

if organism != "*" and organism not in feature_data:
    print("Unknown organism " + organism + " check region data")
    exit(2)

sequence_organisms = {v["sequence"]: k for _, (k, v) in enumerate(feature_data.items())}
organism_sequences = {k: v["sequence"] for _, (k, v) in enumerate(feature_data.items())}

for ok in feature_data:
    org = feature_data[ok]
    for gk in org["genes"]:
        gene = org["genes"][gk]
        gene["region_data"] = [0 for _ in range(gene["start"], gene["stop"] + 1)]

# sequence_id = feature_data[organism]["sequence"]
# genes = feature_data[organism]["genes"]
# region_data = {k: [0 for _ in range(v["start"], v["stop"] + 1)] for i, (k, v) in enumerate(genes.items())}

with gzip.open(pileup, "rt") as f:
    while line := f.readline().rstrip():
        if line.startswith("#") == False:
            elements = line.split("\t")

            if organism == "*" or organism == sequence_organisms[elements[0]]:
                pos = int(elements[1])
                infos = elements[7].split(";")
                for i in infos:
                    v = i.split("=")
                    if i.startswith("DP="):
                        genes = feature_data[sequence_organisms[elements[0]]]["genes"]
                        for g in genes:
                            gene = genes[g]
                            if genes[g]["start"] <= pos and pos <= genes[g]["stop"]:
                                gene["region_data"][pos - genes[g]["start"]] = int(v[1])

print("sample\torganism\tgene\tmeandepth\tmediandepth\tmindepth\tmaxdepth\tstdev\tbases\tseen\tcoverage")
for ok in feature_data:
    if organism == "*" or organism == ok:
        genes = feature_data[ok]["genes"]
        for gk in genes:
            region_data = genes[gk]["region_data"]

            mean = np.mean(region_data)
            median = np.median(region_data)
            min = np.min(region_data)
            max = np.max(region_data)
            stdev = np.std(region_data)
            read_count = sum([1 for dp in region_data if dp != 0])
            width = genes[gk]["stop"] - genes[gk]["start"] + 1

            print(
                "{}\t{}\t{}\t{:.1f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.2f}\t{:.0f}\t{:.0f}\t{:.2f}".format(
                    sample,
                    ok,
                    gk,
                    mean,
                    median,
                    min,
                    max,
                    stdev,
                    width,
                    read_count,
                    read_count / width * 100.0,
                )
            )
