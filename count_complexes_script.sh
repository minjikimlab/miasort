#!/bin/bash
#SBATCH --job-name=count_complexes_parallel
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=count_complexes_parallel.txt
#SBATCH --mem=100g
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --profile=Task

## DO NOT CHANGE ANYTHING BELOW

my_job_header

#### LOAD MODULES / SOURCE ENVIRONMENT HERE
cd /nfs/turbo/umms-minjilab/mia-sort/


######  YOUR PROGRAM CODE BELOW
/nfs/turbo/umms-minjilab/mia-sort/env/bin/python count_complexes.py