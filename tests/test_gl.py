#!/usr/bin/env python3

import miasort

miasort.abc_sort("../test_input.complexes",
                 "../test_input.domains",
                 "AtoC;CtoA;AandC;Bcentered;BtoA;BtoC",
                 out_dir="../20240819_test/test_folder_syn_6000_20240826",
                 anchor_option="yes_complete")

miasort.multiple_sort("../test_input.complexes",
                      "../test_input_abc.domains",
                      out_dir="../20240819_test/test_folder_syn_AandBandC_6000_20240826",
                      anchor_option="yes_complete",
                      subplots_margins=(0.65, 0.15, 0.9))

miasort.unlimited_multiple_sort("../test_input.complexes",
                                "chr3:100000-108000;chr3:150000-155000;chr3:300000-308000;chr3:420000-428000",
                                "yes;no;yes;yes",
                                out_dir="../20240819_test/test_folder_syn_multiple_6000_20240826",
                                anchor_option="yes_complete",
                                subplots_margins=(0.65, 0.15, 0.9))
