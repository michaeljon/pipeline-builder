#!/usr/bin/env python

import sys
import json
import gzip
from os.path import exists, expandvars


if len(sys.argv) != 5:
    print("usage: " + sys.argv[0] + " <sample> <organism> <depth> <destpath>")
    exit(1)

feature_data = {}
with open("hcov-regions.json") as f:
    feature_data = json.load(f)

sample = sys.argv[1]
organism = sys.argv[2]
depth_file = sys.argv[3]
dest = sys.argv[4]

if organism != "*" and organism not in feature_data:
    print("Unknown organism " + organism + " check region data")
    exit(2)

if exists(expandvars(depth_file)) == False:
    print("Unable to location depth file")
    exit(3)

if exists(expandvars(dest)) == False:
    print("Unable to location destination")
    exit(4)

sequence_organisms = {v["sequence"]: k for _, (k, v) in enumerate(feature_data.items())}
organism_sequences = {k: v["sequence"] for _, (k, v) in enumerate(feature_data.items())}

for ok in feature_data:
    org = feature_data[ok]
    org["region_data"] = [0 for _ in range(org["region"]["start"], org["region"]["stop"] + 1)]
    for gk in org["genes"]:
        gene = org["genes"][gk]
        gene["region_data"] = [0 for _ in range(gene["start"], gene["stop"] + 1)]

with gzip.open(depth_file, "rt") as f:
    while line := f.readline().rstrip():
        if line.startswith("#") == False:
            elements = line.split("\t")

            if organism == "*" or organism == sequence_organisms[elements[0]]:
                pos = int(elements[1])

                target = feature_data[sequence_organisms[elements[0]]]
                depth = int(elements[2])

                target["region_data"][pos - 1] = depth

                genes = target["genes"]
                for g in genes:
                    gene = genes[g]
                    if genes[g]["start"] <= pos and pos <= genes[g]["stop"]:
                        gene["region_data"][pos - genes[g]["start"]] = depth

to_process = [k for k in feature_data if organism == "*" or organism == k]

for ok in to_process:
    region_data = feature_data[ok]["region_data"]
    with open(dest + "/" + sample + "_" + ok + "__all" + ".tsv", "w") as f:
        f.write("sample\torganism\tgene\tposition\tdepth\n")
        for position in range(len(region_data)):
            target = feature_data[ok]
            genes = target["genes"]

            gene_name = "<*>"
            for g in genes:
                gene = genes[g]
                if genes[g]["start"] <= position + 1 and position + 1 <= genes[g]["stop"]:
                    gene_name = g

            f.write("{}\t{}\t{}\t{:.0f}\t{:.0f}\n".format(sample, ok, gene_name, position + 1, region_data[position]))

    genes = feature_data[ok]["genes"]
    for gk in genes:
        with open(dest + "/" + sample + "_" + ok + "_" + gk + ".tsv", "w") as f:
            f.write("sample\torganism\tgene\tpos_in_gene\tpos_in_genome\tdepth\n")

            region_data = genes[gk]["region_data"]
            start_pos = genes[gk]["start"]
            for position in range(len(region_data)):
                f.write(
                    "{}\t{}\t{}\t{:.0f}\t{:.0f}\t{:.0f}\n".format(
                        sample, ok, gk, position + 1, start_pos + position, region_data[position]
                    )
                )
