#!/usr/bin/env python

import sys
import os

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from paired_end import *

panel_choices = [
    "MN908947.3",
    "AF304460.1",
    "AY597011.2",
    "AY567487.2",
    "AY585228.1",
    "NC_038311.1",
    "NC_038312.1",
    "NC_038878.1",
    "panel",
]

panel_choice_help = (
    "'Chromosome' name from reference assembly "
    + "(MN908947.3, sars-cov-2), "
    + "(AF304460.1, hcov-229e), "
    + "(AY597011.2, hcov-hku1), "
    + "(AY567487.2, hcov-nl63), "
    + "(AY585228.1, hcov-oc43), "
    + "(NC_038311.1, hrv-a), "
    + "(NC_038312.1, hrv-b), "
    + "(NC_038878.1, hrv-c), "
    + "(panel, combined panel of all organisms)"
)


if __name__ == "__main__":
    main(panel_choices, panel_choice_help)
