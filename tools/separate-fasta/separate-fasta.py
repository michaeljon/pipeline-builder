#!/usr/bin/env python

import argparse

import re
from argparse import Namespace
from os import mkdir, path
from os.path import basename, exists, expandvars, join, splitext
from typing import Any, Dict, List, Sequence, Tuple

# This script reads a FASTA file containing data for multiple organisms and writes
# the data for each organism out as its own FASTA file.

organism_names = {
    # HCoV
    "AF304460.1": "hcov-229e",
    "AY597011.2": "hcov-hku1",
    "AY567487.2": "hcov-nl63",
    "AY585228.1": "hcov-oc43",
    # SARS-CoV-2
    "MN908947.3": "sars-cov-2",
    # HRV
    "NC_001617.1": "hrv-a",
    "NC_038311.1": "hrv-a1",
    "NC_001490.1": "hrv-b14",
    "NC_038312.1": "hrv-b3",
    "NC_009996.1": "hrv-c",
}


def define_arguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        action="store",
        metavar="INPUT_FILE",
        dest="input_file",
        help="Path to FASTA file to separate",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        action="store",
        metavar="OUTPUT_DIR",
        dest="output_dir",
        help="Directory to write output FASTAs to",
    )

    return parser.parse_args()


def expand_paths(opts: Namespace) -> Dict[str, Any]:
    options = vars(opts)

    for opt in ["input_file", "output_dir"]:
        options[opt] = expandvars(options[opt])

    return options


def verify_options(options: Dict[str, Any]):
    if exists(options["input_file"]) == False:
        print("Input file not found at {PATH}".format(PATH=options["input_file"]))
        quit(1)

    if exists(options["output_dir"]) == False:
        print("Output directory not found at {PATH}, creating...".format(PATH=options["output_dir"]))
        mkdir(options["output_dir"])


# Reads a FASTA file and returns a sequence of (organism_id, section_contants) pairs.
def read_fasta_sections(fasta_path: str) -> Sequence[Tuple[str, str]]:
    with open(fasta_path, "r") as fasta:
        sections: List[Tuple[str, str]] = []

        current_organism_id = ""
        current_section_contents = ""

        for line in fasta.readlines():
            if line.startswith(">"):
                if current_section_contents != "":
                    sections.append((current_organism_id, current_section_contents))
                    current_section_contents = ""

                for k in organism_names.keys():
                    if k in line:
                        current_organism_id = k

                line = (
                    ">"
                    + path.basename(fasta_path).replace(".consensus.fa", "")
                    + " | "
                    + current_organism_id
                    + " ("
                    + organism_names[current_organism_id]
                    + ")\n"
                )

            current_section_contents += line

        sections.append((current_organism_id, current_section_contents))

        return sections


# Takes a file name and returns (file_stem, [extensions]).
# Ex. "my_data.consensus.fa" => ("my_data", [".consensus", ".fa"])
def separate_extensions(file_name: str) -> Tuple[str, Sequence[str]]:
    file_stem = file_name
    extensions = []
    while True:
        (new_stem, ext) = splitext(file_stem)
        if new_stem == file_stem:
            extensions.reverse()
            return (file_stem, extensions)
        else:
            file_stem = new_stem
            extensions.append(ext)


if __name__ == "__main__":
    opts = define_arguments()
    options = expand_paths(opts)
    verify_options(options)

    (file_prefix, extensions) = separate_extensions(basename(options["input_file"]))

    for organism_id, section_contents in read_fasta_sections(options["input_file"]):
        file_name = "{PREFIX}.{ORGANISM_NAME}.{EXTENSIONS}".format(
            PREFIX=file_prefix,
            ORGANISM_NAME=organism_names[organism_id],
            EXTENSIONS=".".join(["consensus", "fa"]),
        )
        output_path = join(options["output_dir"], file_name)

        with open(output_path, "w") as output_file:
            print(f"Writing {output_path}")
            output_file.write(section_contents)
