#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

# Record the start time
start_time=$(date +%s)

time cargo run --release --example bed-intersect -- \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno test-july-2.bedte > coitrees_output.bed

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running coitrees intersect: $duration seconds"