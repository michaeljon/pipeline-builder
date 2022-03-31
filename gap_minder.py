#!/usr/bin/env python

import csv
import argparse

from argparse import Namespace
from typing import Dict, Any, List
from os.path import exists, expandvars

OptionsDict = Dict[str, Any]
FieldsDict = Dict[str, Any]
CoverageData = list[int]


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
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
        "-f",
        "--feature-file",
        required=True,
        action="store",
        metavar="FEATURE-FILE",
        dest="featureFile",
        help="Full path to genomic feature file",
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
        default=500,
        help="Minimum depth to trigger gap detection",
    )

    return parser.parse_args()


def fixupPathOptions(opts: Namespace) -> OptionsDict:
    options = vars(opts)

    if options["outputFile"] == "-":
        options["outputFile"] = "/dev/stdout"

    for opt in ["coverageFile", "featureFile", "outputFile"]:
        options[opt] = expandvars(options[opt])

    return options


def verifyOptions(options: OptionsDict):
    for opt in ["coverageFile", "featureFile"]:
        if exists(options[opt]) == False:
            print("Unable to find your --{OPTION} file at {PATH}".format(OPTION=opt, PATH=options[opt]))
            quit(1)


def fieldToDict(type: str, field: str) -> FieldsDict:
    cells = {cell.split("=")[0]: cell.split("=")[1] for cell in [cell for cell in field.split(";")]}

    if "gene" not in cells.keys():
        cells["gene"] = type

    if "ID" not in cells.keys():
        cells["ID"] = ""

    return cells


def loadFeatureData(featureFile: str):
    featureData = []

    with open(featureFile, "r") as f:
        csvFile = csv.reader([row for row in f if not row.startswith("#")], delimiter="\t")

        for line in csvFile:
            if line[2] != "CDS":
                featureData.append(
                    {
                        "type": line[2],
                        "start": int(line[3]),
                        "end": int(line[4]),
                        "length": int(line[4]) - int(line[3]) + 1,
                        "info": fieldToDict(line[2], line[8]),
                    }
                )

    return featureData


def loadCoverageData(coverageFile: str) -> CoverageData:
    coverageData = []
    coverageData.append(-1)

    with open(coverageFile, "r") as f:
        reader = csv.reader([row for row in f if not row.startswith("#")], delimiter="\t")

        for line in reader:
            coverageData.append(int(line[2]))

    return coverageData


def calculateGaps(coverageData: CoverageData, minDepth: int) -> List:
    gaps = []

    lociCount = len(coverageData)
    locus = 1

    while locus < lociCount:
        if coverageData[locus] < minDepth:
            start = locus
            basesLost = 0

            while locus < lociCount and coverageData[locus] < minDepth:
                basesLost += minDepth - coverageData[locus]
                locus += 1

            end = locus - 1
            gaps.append({"start": start, "end": end, "averageDelta": int(basesLost / (end - start + 1))})
        else:
            locus += 1

    return gaps


def overlay(gapData: List, featureData: List) -> List:
    partial = (
        lambda gap, feature: gap["start"] >= feature["start"]
        and gap["start"] <= feature["end"]
        and gap["end"] > feature["end"]
    )
    complete = lambda gap, feature: (gap["start"] >= feature["start"] and gap["start"] <= feature["end"]) and (
        gap["end"] >= feature["start"] and gap["end"] <= feature["end"]
    )

    overlays = []

    for gap in gapData:
        partialOverlays = list(filter(lambda feature: partial(gap, feature), featureData))
        completeOverlays = list(filter(lambda feature: complete(gap, feature), featureData))

        if len(partialOverlays) > 0 or len(completeOverlays) > 0:
            overlays.append(
                {
                    "gapStart": gap["start"],
                    "gapEnd": gap["end"],
                    "gapSize": gap["end"] - gap["start"] + 1,
                    "averageDelta": gap["averageDelta"],
                    "partialOverlays": [
                        {
                            "gene": partial["info"]["gene"],
                            "id": partial["info"]["ID"],
                            "type": partial["type"],
                            "geneStart": partial["start"],
                            "geneEnd": partial["end"],
                        }
                        for partial in partialOverlays
                    ],
                    "completeOverlays": [
                        {
                            "gene": complete["info"]["gene"],
                            "id": complete["info"]["ID"],
                            "type": complete["type"],
                            "geneStart": complete["start"],
                            "geneEnd": complete["end"],
                        }
                        for complete in completeOverlays
                    ],
                }
            )

    return overlays


def writeOutput(overlays: List, coverageFile: str):
    uniqueOverlayKeys = set()
    uniqueOverlays = list()

    for overlay in overlays:
        key = str(overlay["gapStart"]) + "_" + str(overlay["gapEnd"])

        if key not in uniqueOverlayKeys:
            uniqueOverlayKeys.add(key)
            uniqueOverlays.append(overlay)

    with open(coverageFile, "w") as f:
        writer = csv.DictWriter(
            f, ["gapStart", "gapEnd", "gapSize", "averageDelta", "geneStart", "geneEnd", "gene", "id", "overlapType"]
        )
        writer.writeheader()

        for gap in uniqueOverlays:
            for partial in gap["partialOverlays"]:
                writer.writerow(
                    {
                        "gapStart": gap["gapStart"],
                        "gapEnd": gap["gapEnd"],
                        "gapSize": gap["gapSize"],
                        "averageDelta": gap["averageDelta"],
                        "geneStart": partial["geneStart"],
                        "geneEnd": partial["geneEnd"],
                        "gene": partial["gene"],
                        "id": partial["id"],
                        "overlapType": "partial",
                    }
                )

            for complete in gap["completeOverlays"]:
                writer.writerow(
                    {
                        "gapStart": gap["gapStart"],
                        "gapEnd": gap["gapEnd"],
                        "gapSize": gap["gapSize"],
                        "averageDelta": gap["averageDelta"],
                        "geneStart": complete["geneStart"],
                        "geneEnd": complete["geneEnd"],
                        "gene": complete["gene"],
                        "id": complete["id"],
                        "overlapType": "complete",
                    }
                )


def main():
    options = fixupPathOptions(defineArguments())
    verifyOptions(options)

    # load the feature data
    featureData = loadFeatureData(options["featureFile"])

    # load the coverage file
    coverageData = loadCoverageData(options["coverageFile"])

    # locate the coverage gaps based on minimum depth
    gapData = calculateGaps(coverageData, int(options["minDepth"]))

    # overlay the gaps on the feature data
    overlays = overlay(gapData, featureData)

    # write the output
    writeOutput(overlays, options["outputFile"])


if __name__ == "__main__":
    main()
