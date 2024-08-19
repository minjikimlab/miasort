#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++11 -o ./bin/sprite2regions ./bin/sprite2regions.cpp -lz

# Record the start time
start_time=$(date +%s)

time ./bin/sprite2regions . ./data/GM12878-SPRITE.byChromosome.clusters ./data/hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running sprite2regions: $duration seconds"
