#!/usr/bin/env python

import sys

contig = "NC_045512.2"
coverage_data = {}

for filename in sys.argv[1:]:
    sample = filename.split('/')[0]
    coverage_data[sample] = [-1 for cov in range(29904)]

    with open(filename, "r") as sample_coverage:
        for read_depth_data in sample_coverage.readlines():
            read = read_depth_data.split("\t")

            position = int(read[1]) - 1
            coverage = int(read[2])

            coverage_data[sample][position] = coverage

samples = sorted(coverage_data.keys())

with open("coverage.tsv", "w+") as coverage:
    coverage.write("contig\tposition")
    for sample in samples:
        coverage.write("\t" + sample)
    coverage.write("\n")

    for position in range(29903):
        coverage.write(contig + "\t" + str(position + 1))
        for sample in samples:
            coverage.write("\t" + str(coverage_data[sample][position]))
        coverage.write("\n")