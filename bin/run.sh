#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

python ./mia-sort/main.py --path1 ./data/test_input.region --path2 ./data/test_input.domains \
--type abc --graphs AtoC\;CtoA\;AandC\;Bcentered\;BtoA\;BtoC --numfrag_min 2  --anchor_options yes_complete --out_dir ./test_folder_syn_6000

python ./mia-sort/main.py --path1 ./data/test_input.region \
--type multiple --graphs none --numfrag_min 2  \
--region chr3:100000-108000\;chr3:150000-155000\;chr3:300000-308000\;chr3:420000-428000 \
--operation yes\;no\;yes\;yes \
--anchor_options yes_complete --out_dir ./test_folder_syn_multiple_6000