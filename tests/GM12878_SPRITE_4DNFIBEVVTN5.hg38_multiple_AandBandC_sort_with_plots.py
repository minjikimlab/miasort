#!/usr/bin/env python3

import miasort

miasort.multiple_sort("../GM12878_SPRITE_4DNFIBEVVTN5.hg38.ext250bp.g8000bp.region",
                      "../../mia-sort/data/sprite_abc.domains",
                      out_dir="../20240819_test/GM12878_SPRITE_4DNFIBEVVTN5.hg38_multiple_AandBandC_sort_with_plots",
                      extension="natural",
                      anchor_option="yes_top")