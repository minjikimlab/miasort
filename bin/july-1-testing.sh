#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test-july-1.bedte \
--type abc --graphs BtoA\;BtoC\;AtoC\;CtoA\;AandC\;Bcentered --numfrag 2
