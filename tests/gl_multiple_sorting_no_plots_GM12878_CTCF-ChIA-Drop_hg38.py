#!/usr/bin/env python3

import miasort

miasort.multiple_sort("../../GM12878_CTCF-ChIA-Drop_GSE158897.hg38.complexes",
                 "../../GM12878-sprite-three-way-regions.bedte",
                 plot=False,
                 colors="red;#0000FF;#525252",
                 out_dir="../M12878_CTCF-ChIA-Drop-multiple",
                 anchor_option="yes_top")