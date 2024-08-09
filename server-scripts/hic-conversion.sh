#!/bin/bash

#SBATCH --job-name=bedtools_gm_400_line
#SBATCH --mail-user=zhangzzc@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logs/bedtools_gm_400.log

g++ -std=c++11 -o pairs2regions pairs2regions.cpp -lz
./pairs2regions . ./4DNFIACOTIGL.pairs.gz ./hg38.chrom.sizes
