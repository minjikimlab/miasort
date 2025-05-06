#!/usr/bin/env python3

import sys
import os

# Ensure the local miasort module is used
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + "/.."))

import miasort

miasort.abc_sort("./data/test_input.region",
                 "./data/test_input.domains",
                 "AtoC;CtoA;AandC;Bcentered;BtoA;BtoC",
                 histogram=True, # Enable histogram plotting
                 out_dir="./test_folder_syn_6000",
                 anchor_option="yes_complete")