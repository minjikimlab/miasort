#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++11 -o pairs2regions pairs2regions.cpp -lz

# Record the start time
start_time=$(date +%s)

time ./pairs2regions . ./4DNFIACOTIGL.pairs.gz ./hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running pairs2regions: $duration seconds"
