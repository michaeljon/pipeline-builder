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


def loadFeatureData(regionFile: str):
    feature_data = {}

    with open(regionFile) as f:
        feature_data = json.load(f)

    for k in feature_data.keys():
        organism_to_sequence[k] = feature_data[k]["sequence"]
        sequence_to_organism[feature_data[k]["sequence"]] = k

    return feature_data


def loadCoverageData(coverageFile: str, featureData) -> CoverageData:
    coverageData = {}

    for org in featureData.keys():
        start = int(featureData[org]["region"]["start"])
        stop = int(featureData[org]["region"]["stop"])

        # add an extra 1000 bases on to the organism in case of a large insert
        coverageData[org] = [0 for _ in range(start, stop + 1000)]

    with gzip.open(coverageFile, "rt") as f:
        while line := f.readline().rstrip():
            if line.startswith("#") == False:
                elements = line.split("\t")

                organism = sequence_to_organism[elements[0]]
                position = int(elements[1])
                depth = int(elements[2])

                coverageData[organism][position] = depth

    return coverageData


def calculateGaps(coverageData: CoverageData, featureData, minDepth: int) -> GapList:
    gaps = {}

    for org in featureData.keys():
        gaps[org] = []

        lociCount = len(coverageData[org])
        locus = 1

        while locus < lociCount:
            if coverageData[org][locus] < minDepth:
                start = locus
                deltas = []
                depths = []

                while locus < lociCount and coverageData[org][locus] < minDepth:
                    deltas.append(minDepth - coverageData[org][locus])
                    depths.append(coverageData[org][locus])

                    locus += 1

                gaps[org].append(
                    {
                        "start": start,
                        "end": locus - 1,
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
                    }
                )
            else:
                locus += 1

    return gaps


def gatherOverlays(overlays):
    return [
        {
            "gene": overlay["gene"],
            "type": overlay["type"],
            "geneStart": overlay["start"],
            "geneEnd": overlay["end"],
        }
        for overlay in overlays
    ]


def overlay(gapData: GapList, featureData: Any) -> OverlayList:
    overlays = {}

    for org in featureData.keys():
        overlays[org] = []

        partial = (
            lambda gap, feature: gap["start"] >= feature["start"]
            and gap["start"] <= feature["end"]
            and gap["end"] > feature["end"]
        )
        complete = lambda gap, feature: (gap["start"] >= feature["start"] and gap["start"] <= feature["end"]) and (
            gap["end"] >= feature["start"] and gap["end"] <= feature["end"]
        )

        for gap in gapData[org]:
            features = [
                {
                    "start": featureData[org]["genes"][f]["start"],
                    "end": featureData[org]["genes"][f]["stop"],
                    "gene": f,
                    "type": "gene",
                }
                for f in featureData[org]["genes"].keys()
            ]

            partialOverlays = list(filter(lambda feature: partial(gap, feature), features))
            completeOverlays = list(filter(lambda feature: complete(gap, feature), features))

            if len(partialOverlays) > 0 or len(completeOverlays) > 0:
                overlays[org].append(
                    {
                        "organism": org,
                        "gapStart": gap["start"],
                        "gapEnd": gap["end"],
                        "gapSize": gap["end"] - gap["start"] + 1,
                        "averageDepth": gap["averageDepth"],
                        "medianDepth": gap["medianDepth"],
                        "minDepth": gap["minDepth"],
                        "maxDepth": gap["maxDepth"],
                        "stdevDepth": gap["stdevDepth"],
                        "quantileDepth": gap["quantileDepth"],
                        "averageDelta": gap["averageDelta"],
                        "medianDelta": gap["medianDelta"],
                        "minDelta": gap["minDelta"],
                        "maxDelta": gap["maxDelta"],
                        "stdevDelta": gap["stdevDelta"],
                        "quantileDelta": gap["quantileDelta"],
                        "partialOverlays": gatherOverlays(partialOverlays),
                        "completeOverlays": gatherOverlays(completeOverlays),
                    }
                )

    return overlays


def writeRow(writer: csv.DictWriter, sample, gap, overlay, type):
    writer.writerow(
        {
            "sample": sample,
            "organism": gap["organism"],
            "gapStart": gap["gapStart"],
            "gapEnd": gap["gapEnd"],
            "gapSize": gap["gapSize"],
            "averageDepth": gap["averageDepth"],
            "medianDepth": gap["medianDepth"],
            "minDepth": gap["minDepth"],
            "maxDepth": gap["maxDepth"],
            "stdevDepth": gap["stdevDepth"],
            "firstQuartileDepth": gap["quantileDepth"][0],
            "thirdQuartileDepth": gap["quantileDepth"][2],
            "averageDelta": gap["averageDelta"],
            "medianDelta": gap["medianDelta"],
            "minDelta": gap["minDelta"],
            "maxDelta": gap["maxDelta"],
            "stdevDelta": gap["stdevDelta"],
            "firstQuartileDelta": gap["quantileDelta"][0],
            "thirdQuartileDelta": gap["quantileDelta"][2],
            "geneStart": overlay["geneStart"],
            "geneEnd": overlay["geneEnd"],
            "gene": overlay["gene"],
            "overlapType": type,
        }
    )


def writeOutput(overlays: OverlayList, organism: str, coverageFile: str, featureData: Any, sample: str):
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
                "overlapType",
            ],
        )
        writer.writeheader()

        to_process = [k for k in featureData if organism == "*" or organism == k]
        for org in to_process:
            uniqueOverlayKeys = set()
            uniqueOverlays = list()

            for overlay in overlays[org]:
                key = str(overlay["gapStart"]) + "_" + str(overlay["gapEnd"])

                if key not in uniqueOverlayKeys:
                    uniqueOverlayKeys.add(key)
                    uniqueOverlays.append(overlay)

            for gap in uniqueOverlays:
                for overlay in gap["partialOverlays"]:
                    writeRow(writer, sample, gap, overlay, "partial")

                for overlay in gap["completeOverlays"]:
                    writeRow(writer, sample, gap, overlay, "complete")


def main():
    options = fixupPathOptions(defineArguments())
    verifyOptions(options)

    # load the feature data
    featureData = loadFeatureData(options["regionFile"])

    # load the coverage file
    coverageData = loadCoverageData(options["coverageFile"], featureData)

    # locate the coverage gaps based on minimum depth
    gapData = calculateGaps(coverageData, featureData, int(options["minDepth"]))

    # overlay the gaps on the feature data
    overlays = overlay(gapData, featureData)

    # write the output
    writeOutput(overlays, options["organism"], options["outputFile"], featureData, options["sample"])


if __name__ == "__main__":
    main()
