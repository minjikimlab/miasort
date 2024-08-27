#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../GM12878_insitu-Hi-C_4DNFIYECESRC.complexes",
                 "../../GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="../GM12878_insitu-Hi-C_4DNFIYECESRC_conv-loops_stripes",
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")