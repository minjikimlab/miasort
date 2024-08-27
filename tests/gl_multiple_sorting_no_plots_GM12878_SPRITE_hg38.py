#!/usr/bin/env python3

import miasort

miasort.multiple_sort("../../GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes",
                 "../../GM12878-sprite-three-way-regions.bedte",
                 plot=False,
                 colors="red;#FF0000;#525252",
                 out_dir="../GM12878_SPRITE-multiple_no_plots",
                 anchor_option="yes_top")