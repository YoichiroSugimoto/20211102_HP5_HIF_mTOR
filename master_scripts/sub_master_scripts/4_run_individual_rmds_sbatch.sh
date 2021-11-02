#!/bin/bash
#SBATCH --job-name=4_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=3-00:00:0
#SBATCH --mem=140G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch --dependency=afterok:24925462 /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/4_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate run_dpi_v2

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/R/s4-tss-definition

Rscript ../run_rmd.R s4-1-1-tss-definition_v3.rmd
Rscript ../run_rmd.R s4-1-2-tss-definition_v2.rmd

conda activate five_prime_seq_for_VHL_loss_v0.2.1

Rscript ../run_rmd.R s4-2-1-tx-tss-assignment-for-RCC4-v5.rmd
Rscript ../run_rmd.R s4-3-1-extract-tx-information-for-RCC4-v2.rmd

Rscript ../run_rmd.R s4-2-2-tx-tss-assignment-for-786O-v5.rmd
Rscript ../run_rmd.R s4-3-2-extract-tx-information-for-786O-v2.rmd


## Template for nohup command
## nohup Rscript ../run_rmd.R x.rmd &> ../nohup_out/x.nohup.out
##
## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
