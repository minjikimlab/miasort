#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++11 -o ../mia-sort/bin/pairs2regions ../mia-sort/bin/pairs2regions.cpp -lz -lpthread

# Record the start time
start_time=$(date +%s)

time ./bin/pairs2regions . ./data/LHG0035N_0035V_0045V.bsorted.pairs.gz ./data/hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running pairs2regions: $duration seconds"
