#!/usr/bin/env python3

import miasort

miasort.abc_sort("../../GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes",
                 "../../GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="../GM12878_SPRITE_4DNFIBEVVTN5-Drop_conv-loops_stripes",
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")