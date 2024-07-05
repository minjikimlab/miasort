#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 test-july-2.bedte \
--type abc --graphs BtoA\;BtoC\;AtoC\;CtoA\;AandC\;Bcentered --numfrag_min 2  --numfrag_max 2 --out_dir test_folder_4_numfrag