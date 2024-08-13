#!/bin/bash

#SBATCH --job-name=bedtools_hic_1_line
#SBATCH --mail-user=zhangzzc@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logs/bedtools_hic_1.log
#SBATCH --mem=50000m
#SBATCH --time=01-02:00:00

python main.py --path1 4DNFIACOTIGL.pairs.gz.ext250bp.g8000bp.region --path2 anchors-1.bedte \
--type abc --graphs BtoA\;BtoC\;AtoB\;AtoC\;CtoA\;CtoB\;AandC\;Bcentered --numfrag_min 2  --anchor_options no --out_dir test_greatlakes_bedtools_hic_1
