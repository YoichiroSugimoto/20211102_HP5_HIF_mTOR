#!/bin/bash
#SBATCH --job-name=7_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch --dependency=afterok:24925486 /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/7_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/R/s7-HIF-binding-site-definition

Rscript ../run_rmd.R s7-1-curate-ENCODE-pipeline-processed-ChIP-Seq-data.rmd
Rscript ../run_rmd.R s7-2-HIF1A-and-HIF2A-binding-ratio.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
