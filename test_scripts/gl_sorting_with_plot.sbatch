#!/bin/bash

#SBATCH --job-name=gl_sorting_with_plot
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=10
#SBATCH --mem=60g
#SBATCH --gpus=0
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --profile=Task
#SBATCH --out=jet_with_plot.out

## DO NOT CHANGE ANYTHING BELOW

my_job_header

#### LOAD MODULES / SOURCE ENVIRONMENT HERE
cd /nfs/turbo/umms-minjilab/mia-sort/
pwd
conda activate miasort

######  YOUR PROGRAM CODE BELOW
./tests/gl_jet_sorting_with_plots_DP-thymocyte_mm10.py