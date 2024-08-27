#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes",
                 "../../GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="../GM12878_CTCF-ChIA-PET-Drop_conv-loops_stripes",
                 colors="red;#0000FF;#525252",
                 anchor_option="yes_top")