#!/bin/bash

#SBATCH --job-name=bedtools_sprite_100_line
#SBATCH --mail-user=zhangzzc@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./logs/bedtools_sprite_100.log
#SBATCH --mem=12000m
#SBATCH --time=01-02:00:00

python main.py --path1 GM12878-SPRITE.byChromosome.clusters --path2 anchors-100.bedte \
--type abc --graphs BtoA\;BtoC\;AtoB\;AtoC\;CtoA\;CtoB\;AandC\;Bcentered --numfrag_min 2  --anchor_options no --out_dir test_greatlakes_bedtools_sprite_100
