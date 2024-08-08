#!/bin/bash

# similar setting as the GreatLakes testing

# Stop on errors
set -Eeuo pipefail
set -x

# Record the start time
start_time=$(date +%s)

# Run the Python script with time measurement
time python main.py --path1 4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region --path2 anchors-400.bedte \
--type abc --graphs BtoA\;BtoC\;AtoB\;AtoC\;CtoA\;CtoB\;AandC\;Bcentered --numfrag_min 2  --anchor_options no --out_dir test_folder_400_hic

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running ComplexSorter: $duration seconds"
