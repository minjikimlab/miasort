#!/bin/bash

# test runtime on a small BED file

# Stop on errors
set -Eeuo pipefail
set -x

# Record the start time
start_time=$(date +%s)

# sort -k1,1 -k2,2n GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno > in.sorted.bed

# Run the Python script with time measurement
time python ./mia-sort/main.py --path1 ./data/GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 ./data/anchors-5.bedte \
--type abc --graphs AtoC\;CtoA\;AandC\;Bcentered\;BtoA\;BtoC --numfrag_min 2  --anchor_options yes_top --out_dir ./test_folder_unit

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Total time for running ComplexSorter: $duration seconds"
