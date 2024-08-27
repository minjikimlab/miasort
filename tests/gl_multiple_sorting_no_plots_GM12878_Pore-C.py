#!/usr/bin/env python3

import miasort

miasort.multiple_sort("../../GM12878_Pore-C_GSM4490689.hg38.complexes",
                 "../../GM12878-sprite-three-way-regions.bedte",
                 plot=False,
                 colors="red;#FF0000;#525252",
                 out_dir="../GM12878_Pore-C-multiple_no_plots",
                 anchor_option="yes_top")