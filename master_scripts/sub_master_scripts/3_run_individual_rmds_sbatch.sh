#!/bin/bash
#SBATCH --job-name=3_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:0
#SBATCH --mem=56G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/3_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.1.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/s3-alignment-statistics

Rscript ../run_rmd.R s3-1-alignment-statistics.rmd
Rscript ../run_rmd.R s3-2-evaluation-of-deduplication-with-UMI.rmd

conda activate five_prime_seq_for_VHL_loss_v0.2.1

Rscript ../run_rmd.R s3-3-PCA.rmd
Rscript ../run_rmd.R s3-4-ERCC-evaluation.rmd

## Template for nohup command
## nohup Rscript ../run_rmd.R x.rmd &> ../nohup_out/x.nohup.out
##
## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
