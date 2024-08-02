#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-1.bedte > out-coitrees-1.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-5.bedte > out-coitrees-5.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-10.bedte > out-coitrees-10.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-100.bedte > out-coitrees-100.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-200.bedte > out-coitrees-200.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno anchors-400.bedte > out-coitrees-400.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-1.bedte > out-coitrees-1.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-5.bedte > out-coitrees-5.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-10.bedte > out-coitrees-10.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-100.bedte > out-coitrees-100.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-200.bedte > out-coitrees-200.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"

start_time=$(date +%s)
time cargo run --release --example bed-intersect -- -i \
4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region anchors-400.bedte > out-coitrees-400.bed
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Total time: $duration secs"
