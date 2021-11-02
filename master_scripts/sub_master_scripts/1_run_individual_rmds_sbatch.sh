#!/bin/bash
#SBATCH --job-name=1_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --time=3-00:00:0
#SBATCH --mem=192G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/1_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.1.1

## Preprocessing input fastq files
cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/s1-preprocessing-reads

Rscript ../run_rmd.R s1-1-preprocessing_of_fastq-step1.rmd
Rscript ../run_rmd.R s1-2-preprocessing_of_fastq-step2.rmd


## Template for nohup command
## nohup Rscript ../run_rmd.R x.rmd &> ../nohup_out/x.nohup.out
##
## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
