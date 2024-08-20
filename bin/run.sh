#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

python ./mia-sort/main.py --path1 ./data/test_input.region --path2 ./data/test_input.domains \
--type abc --graphs AtoC\;CtoA\;AandC\;Bcentered\;BtoA\;BtoC --numfrag_min 2  --anchor_options yes_complete --out_dir ./test_folder_syn_6000