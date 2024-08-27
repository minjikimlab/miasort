#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../GM12878_Pore-C_GSM4490689.hg38.complexes",
                 "../../GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="../GM12878_Pore-C_GSM4490689-Drop_conv-loops_stripes",
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")