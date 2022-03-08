#!/bin/bash
#SBATCH --job-name=10_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/10_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate hydroxylation_by_JMJD6

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/s10-additional-analysis

Rscript ../run_rmd.R s10-1-validation-of-HP5-by-RT-qPCR.rmd
Rscript ../run_rmd.R s10-2-HP5-dynamic-range.rmd
Rscript ../run_rmd.R s10-3-comparison-with-orthogonal-methods.rmd


## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
