#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../DP-thymocyte_Hi-C_GSE199059.mm10.complexes",
                 "../../DP-thymocyte-jets.bedte",
                 "Bcentered;BtoA;BtoC",
                 plot=False,
                 colors="red;#FF0000;#525252",
                 out_dir="../DP-thymocyte_Hi-C_GSE199059-DP-thymocyte-jets",
                 anchor_option="yes_top")