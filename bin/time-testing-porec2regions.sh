#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -std=c++17 -o ./bin/porec2complexes ./bin/porec2complexes.cpp -lz

# Record the start time
start_time=$(date +%s)

./bin/porec2complexes ./data/ ./data/GSM4490689_GM12878_DpnII_100db_GRCh38_bwa_0.7.17_sensitive_GIABhiconf_whatshap_0.19.c01520_default.fragment_alignments_first_1k.csv ./data/hg38.chrom.sizes \
./data/GSM4490689_GM12878_DpnII_100db_GRCh38_bwa_0.7.17_sensitive_GIABhiconf_whatshap_0.19.c01520_default.fragment_alignments

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running porec2regions: $duration seconds"
