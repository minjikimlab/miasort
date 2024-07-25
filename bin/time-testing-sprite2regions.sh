#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++11 -o sprite2regions sprite2regions.cpp -lz

# Record the start time
start_time=$(date +%s)

time ./sprite2regions . ./GM12878-SPRITE.byChromosome.clusters ./hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running sprite2regions: $duration seconds"
