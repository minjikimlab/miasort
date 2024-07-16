#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

g++ -o pairs2regions pairs2regions.cpp -lz

./pairs2regions . ./LHG0035N_0035V_0045V.bsorted.pairs.gz ./hg38.chrom.sizes
