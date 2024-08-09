#!/bin/bash

#SBATCH --job-name=bedtools_gm_10_line
#SBATCH --mail-user=zhangzzc@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logs/bedtools_gm_10.log

python main.py --path1 GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno --path2 anchors-10.bedte \
--type abc --graphs BtoA\;BtoC\;AtoB\;AtoC\;CtoA\;CtoB\;AandC\;Bcentered --numfrag_min 2  --anchor_options no --out_dir test_greatlakes_bedtools_gm_10
