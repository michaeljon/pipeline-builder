#!/usr/bin/env python

import sys
import os

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from paired_end import *

if __name__ == "__main__":
    main()
