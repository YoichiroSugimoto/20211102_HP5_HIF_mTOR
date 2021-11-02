#!/bin/bash
#SBATCH --job-name=5_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch --dependency=afterok:24925466 /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/5_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/R/s5-tss-validation

Rscript ../run_rmd.R s5-1-validation-of-TSS.rmd

## Template for nohup command
## nohup Rscript ../run_rmd.R x.rmd &> ../nohup_out/x.nohup.out
##
## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
