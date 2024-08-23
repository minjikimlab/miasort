#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++17 -o ./bin/sprite2complexes ./bin/sprite2complexes.cpp -lz

# Record the start time
start_time=$(date +%s)

time ./bin/sprite2complexes . ./data/4DNFIBEVVTN5.clusters.gz ./data/hg38.chrom.sizes

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running sprite2regions: $duration seconds"
