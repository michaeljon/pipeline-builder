#!/usr/bin/env python

import sys
import os

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from paired_end import *

panel_choices = ["NC_038311.1", "NC_038312.1", "NC_038878.1", "panel"]

panel_choice_help = (
    "'Chromosome' name from reference assembly "
    + "(NC_038311.1, hrv-a), "
    + "(NC_038312.1, hrv-b), "
    + "(NC_038878.1, hrv-c), "
    + "(panel, combined panel of all organisms)"
)


if __name__ == "__main__":
    main(panel_choices, panel_choice_help)
