#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../GM12878_CTCF-ChIA-Drop_GSE158897.hg38.complexes",
                 "../../GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="../GM12878_CTCF-ChIA-Drop_conv-loops_stripes",
                 colors="red;#0000FF;#525252",
                 anchor_option="yes_top")