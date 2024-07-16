#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -o pairs2regions pairs2regions.cpp -lz

# Record the start time
start_time=$(date +%s)

time ./pairs2regions . ./LHG0035N_0035V_0045V.bsorted.pairs.gz ./hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running pairs2regions: $duration seconds"
