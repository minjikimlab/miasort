#!/usr/bin/env python3

import miasort

miasort.abc_sort("./data/test_input.region",
                 "./data/test_input.domains",
                 "AtoC;CtoA;AandC;Bcentered;BtoA;BtoC",
                 out_dir="./test_folder_syn_6000",
                 anchor_option="yes_complete")