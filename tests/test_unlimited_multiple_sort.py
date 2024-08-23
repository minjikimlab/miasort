#!/usr/bin/env python3

import miasort

miasort.unlimited_multiple_sort("./data/test_input.region",
                                "chr3:100000-108000;chr3:150000-155000;chr3:300000-308000;chr3:420000-428000",
                                "yes;no;yes;yes",
                                out_dir="./test_folder_syn_multiple_6000",
                                anchor_option="yes_complete",
                                subplots_margins=(0.65, 0.15, 0.9))