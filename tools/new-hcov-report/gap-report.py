#!/usr/bin/env python

import csv
import argparse
from statistics import median, mean, pstdev, quantiles

from argparse import Namespace
import gzip
import sys
import json
from typing import Dict, Any, List
from os.path import exists, expandvars

# This report generates a "gap analysis" for the given sample. It does this by
# reading through the depth.gz (produced by samtools). For each organism we read
# the resulting depths across that's organism's features. Where a given depth is
# reported as lower than the --min-depth a gap is reported.
#
# For each reported gap the organism, gene, and intersection type is reported. The
# intersection describes whether the gap is completely enclosed by the feature (complete),
# or if it intersects multiple features / straddles an edge (partial). Each gap is
# then reported with the following stats:
#
# the gap's location and size: "gapStart", "gapEnd", "gapSize"
# actual read depth across the gap: min, max, avg, median, stdev, 1st and 3rd quartile
# depth difference from min-depth across the gap: min, max, avg, median, stdev, 1st and 3rd quartile
# gene /feature information: "gene", "geneStart", "geneEnd"
# the gap / feature intersection type


OptionsDict = Dict[str, Any]
FieldsDict = Dict[str, Any]
CoverageData = Dict[str, list[int]]
GapList = Dict[str, List]
OverlayList = Dict[str, List]

sequence_to_organism = {}
organism_to_sequence = {}


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--sample",
        required=False,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="Short name for sample",
    )
    parser.add_argument(
        "-c",
        "--coverage-file",
        required=True,
        action="store",
        metavar="COVERAGE-FILE",
        dest="coverageFile",
        help="Full path to file holding coverage data for sample",
    )
    parser.add_argument(
        "-r",
        "--region-file",
        required=True,
        action="store",
        metavar="REGION-FILE",
        dest="regionFile",
        default="hcov-regions.json",
        help="Full path to hcov-regions.json",
    )
    parser.add_argument(
        "-f",
        "--organism",
        required=True,
        action="store",
        metavar="ORGANISM",
        dest="organism",
        default="*",
        help="Name of target organism, or '*' to process all",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        action="store",
        metavar="OUTPUT-FILE",
        dest="outputFile",
        help="Full path to output file or `-` for stdout",
    )
    parser.add_argument(
        "-d",
        "--min-depth",
        required=False,
        action="store",
        metavar="MIN-DEPTH",
        dest="minDepth",
        default=1,
        help="Minimum depth to trigger gap detection",
    )

    return parser.parse_args()


def fixupPathOptions(opts: Namespace) -> OptionsDict:
    options = vars(opts)

    if options["outputFile"] == "-":
        options["outputFile"] = "/dev/stdout"

    for opt in ["coverageFile", "regionFile", "outputFile"]:
        options[opt] = expandvars(options[opt])

    return options


def verifyOptions(options: OptionsDict):
    for opt in ["coverageFile", "regionFile"]:
        if exists(options[opt]) == False:
            print("Unable to find your --{OPTION} file at {PATH}".format(OPTION=opt, PATH=options[opt]))
            quit(1)


def loadFeatureData(regionFile: str, organism: str):
    feature_data = {}

    with open(regionFile) as f:
        feature_data = json.load(f)

    for k in feature_data.keys():
        organism_to_sequence[k] = feature_data[k]["sequence"]
        sequence_to_organism[feature_data[k]["sequence"]] = k

    if organism != "*":
        feature_data = {organism: feature_data[organism]}

    return feature_data


def loadCoverageData(coverageFile: str, organism: str, featureData) -> CoverageData:
    coverageData = {}

    for org in featureData.keys():
        start = int(featureData[org]["region"]["start"])
        stop = int(featureData[org]["region"]["stop"])

        # add an extra 1000 bases on to the organism in case of a large insert
        coverageData[org] = [0 for _ in range(start, stop + 1000)]

    with gzip.open(coverageFile, "rt") as f:
        depthReader = csv.DictReader(f, ["accession", "position", "depth"], delimiter="\t")

        for row in depthReader:
            org = sequence_to_organism[row["accession"]]
            position = int(row["position"])
            depth = int(row["depth"])

            coverageData[org][position] = depth

    return coverageData


def calculateGaps(coverageData: CoverageData, sample, organism_data, minDepth: int) -> GapList:
    gaps = {}

    for organism in organism_data.keys():
        gaps[organism] = []

        for gene, feature in organism_data[organism]["genes"].items():

            locus = int(feature["start"])
            locusEnd = int(feature["stop"])

            start = -1
            deltas = []
            depths = []

            while locus <= locusEnd:
                if coverageData[organism][locus] < minDepth:
                    start = locus

                    while locus <= locusEnd and coverageData[organism][locus] < minDepth:
                        deltas.append(minDepth - coverageData[organism][locus])
                        depths.append(coverageData[organism][locus])

                        locus += 1

                    quantileDepth = quantiles(depths, n=4, method="inclusive") if len(depths) > 2 else [0, 0, 0]
                    quantileDelta = quantiles(deltas, n=4, method="inclusive") if len(deltas) > 2 else [0, 0, 0]

                    gaps[organism].append(
                        {
                            "sample": sample,
                            "organism": organism,
                            "gene": gene,
                            "gapStart": start,
                            "gapEnd": locus - 1,
                            "averageDepth": round(mean(depths), 2),
                            "medianDepth": median(depths),
                            "minDepth": min(depths),
                            "maxDepth": max(depths),
                            "stdevDepth": round(pstdev(depths), 2) if len(depths) > 2 else 0,
                            "firstQuartileDepth": quantileDepth[0],
                            "thirdQuartileDepth": quantileDepth[2],
                            "averageDelta": round(mean(deltas), 2),
                            "medianDelta": median(deltas),
                            "minDelta": min(deltas),
                            "maxDelta": max(deltas),
                            "stdevDelta": round(pstdev(deltas), 2) if len(deltas) > 2 else 0,
                            "firstQuartileDelta": quantileDelta[0],
                            "thirdQuartileDelta": quantileDelta[2],
                            "geneStart": int(feature["start"]),
                            "geneEnd": int(feature["stop"]),
                        }
                    )

                    deltas = []
                    depths = []
                else:
                    locus += 1

            if len(deltas) > 0 or len(depths) > 0:
                # append the last one
                quantileDepth = quantiles(depths, n=4, method="inclusive") if len(depths) > 2 else [0, 0, 0]
                quantileDelta = quantiles(deltas, n=4, method="inclusive") if len(deltas) > 2 else [0, 0, 0]

                gaps[organism].append(
                    {
                        "organism": organism,
                        "gene": gene,
                        "gapStart": start,
                        "gapEnd": locus - 1,
                        "averageDepth": round(mean(depths), 2),
                        "medianDepth": median(depths),
                        "minDepth": min(depths),
                        "maxDepth": max(depths),
                        "stdevDepth": round(pstdev(depths), 2) if len(depths) > 2 else 0,
                        "quantileDepth": quantiles(depths, n=4, method="inclusive") if len(depths) > 2 else [0, 0, 0],
                        "averageDelta": round(mean(deltas), 2),
                        "medianDelta": median(deltas),
                        "minDelta": min(deltas),
                        "maxDelta": max(deltas),
                        "stdevDelta": round(pstdev(deltas), 2) if len(deltas) > 2 else 0,
                        "quantileDelta": quantiles(deltas, n=4, method="inclusive") if len(deltas) > 2 else [0, 0, 0],
                        "geneStart": int(feature["start"]),
                        "geneEnd": int(feature["stop"]),
                    }
                )

    return gaps


# def writeRow(writer: csv.DictWriter, sample, gap, overlay, type):
#     writer.writerow(
#         {
#             "sample": sample,
#             "organism": gap["organism"],
#             "gapStart": gap["gapStart"],
#             "gapEnd": gap["gapEnd"],
#             "gapSize": gap["gapSize"],
#             "averageDepth": gap["averageDepth"],
#             "medianDepth": gap["medianDepth"],
#             "minDepth": gap["minDepth"],
#             "maxDepth": gap["maxDepth"],
#             "stdevDepth": gap["stdevDepth"],
#             "firstQuartileDepth": gap["quantileDepth"][0],
#             "thirdQuartileDepth": gap["quantileDepth"][2],
#             "averageDelta": gap["averageDelta"],
#             "medianDelta": gap["medianDelta"],
#             "minDelta": gap["minDelta"],
#             "maxDelta": gap["maxDelta"],
#             "stdevDelta": gap["stdevDelta"],
#             "firstQuartileDelta": gap["quantileDelta"][0],
#             "thirdQuartileDelta": gap["quantileDelta"][2],
#             "gene": overlay["gene"],
#             "geneStart": overlay["geneStart"],
#             "geneEnd": overlay["geneEnd"],
#         }
#     )


def writeOutput(gap_data: GapList, organism: str, coverageFile: str, organism_data: Any):
    with open(coverageFile, "w") as f:
        writer = csv.DictWriter(
            f,
            [
                "sample",
                "organism",
                "gapStart",
                "gapEnd",
                "gapSize",
                "averageDepth",
                "medianDepth",
                "minDepth",
                "maxDepth",
                "stdevDepth",
                "firstQuartileDepth",
                "thirdQuartileDepth",
                "averageDelta",
                "medianDelta",
                "minDelta",
                "maxDelta",
                "stdevDelta",
                "firstQuartileDelta",
                "thirdQuartileDelta",
                "geneStart",
                "geneEnd",
                "gene",
            ],
        )
        writer.writeheader()

        for organism in organism_data:
            writer.writerows(gap_data[organism])


def main():
    options = fixupPathOptions(defineArguments())
    verifyOptions(options)

    # load the feature data
    organism_data = loadFeatureData(options["regionFile"], options["organism"])

    # load the coverage file
    coverage_data = loadCoverageData(options["coverageFile"], options["organism"], organism_data)

    # locate the coverage gaps based on minimum depth
    gap_data = calculateGaps(coverage_data, options["sample"], organism_data, int(options["minDepth"]))

    # write the output
    writeOutput(gap_data, options["organism"], options["outputFile"], organism_data)


if __name__ == "__main__":
    main()
