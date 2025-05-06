#!/bin/bash
#SBATCH --job-name=cython_memory_no_plot_repeats10_CTCF-ChIA-PET
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/cr527_test_memory/no_plot_repeats10_CTCF-ChIA-PET_hg38.txt
#SBATCH --mem=50g
#SBATCH --time=02:00:00
#SBATCH --profile=Task

## DO NOT CHANGE ANYTHING BELOW

my_job_header

#### LOAD MODULES / SOURCE ENVIRONMENT HERE
cd /nfs/turbo/umms-minjilab/mia-sort/
#source /nfs/turbo/umms-minjilab/mia-sort/env/bin/activate
#conda activate miasort
#export PYTHONPATH=$(pwd)

######  YOUR PROGRAM CODE BELOW
/nfs/turbo/umms-minjilab/mia-sort/env/bin/python -m memory_profiler ./tests/gl_cr527_profiling_index_memory_10repeats_no_plot.py 