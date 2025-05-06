#!/bin/bash
#SBATCH --job-name=test_repeat100_time_profiling
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/cr527_test_time/no_plot_repeats100_CTCF-ChIA-PET_hg38.txt
#SBATCH --mem=100g
#SBATCH --time=02:00:00
#SBATCH --profile=Task

## DO NOT CHANGE ANYTHING BELOW

my_job_header

#### LOAD MODULES / SOURCE ENVIRONMENT HERE
cd /nfs/turbo/umms-minjilab/mia-sort/

######  YOUR PROGRAM CODE BELOW
kernprof -l -v ./tests/gl_cr527_profiling_index_time_100_repeats_no_plot.py 