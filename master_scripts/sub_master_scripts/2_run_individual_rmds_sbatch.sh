#!/bin/bash
#SBATCH --job-name=2_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --time=3-00:00:0
#SBATCH --mem=192G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/2_run_individual_rmds_sbatch.sh
## or
# sbatch --dependency=afterok:24925457 /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/2_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.1.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/R/s2-alignment-of-reads

Rscript ../run_rmd.R s2-1-alignment_of_reads_with_STAR_2-pass-mode.rmd
Rscript ../run_rmd.R s2-2-process-aligned-reads.rmd

## Template for nohup command
## nohup Rscript ../run_rmd.R x.rmd &> ../nohup_out/x.nohup.out
##
## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
