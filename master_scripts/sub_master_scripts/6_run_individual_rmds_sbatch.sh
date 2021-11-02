#!/bin/bash
#SBATCH --job-name=6_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch --dependency=afterok:24925470 /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/6_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/s6-differential-expression-and-tss-usage

Rscript ../run_rmd.R s6-1-gene-level-mRNA-abundance-change-analysis.rmd
Rscript ../run_rmd.R s6-2-tx-level-mRNA-abundance-change-and-diff-TSS.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
