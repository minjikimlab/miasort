#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++11 -o porec2regions porec2regions.cpp -lz

# Record the start time
start_time=$(date +%s)

./porec2regions . ./GSM4490689_GM12878_DpnII_100db_GRCh38_bwa_0.7.17_sensitive_GIABhiconf_whatshap_0.19.c01520_default.fragment_alignments.csv.gz

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running porec2regions: $duration seconds"
