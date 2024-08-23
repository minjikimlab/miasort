#!/usr/bin/env python3

import miasort

miasort.multiple_sort("./data/test_input.region",
                      "./data/test_input_abc.domains",
                      out_dir="./test_folder_syn_AandBandC_6000",
                      anchor_option="yes_complete",
                      subplots_margins=(0.65, 0.15, 0.9))