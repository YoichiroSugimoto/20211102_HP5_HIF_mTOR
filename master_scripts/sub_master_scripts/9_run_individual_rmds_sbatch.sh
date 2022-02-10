#!/bin/bash
#SBATCH --job-name=9_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/9_run_individual_rmds_sbatch.sh
# or
# sbatch --dependency=afterok:25498947 /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/9_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/s9-integrative-analysis

Rscript ../run_rmd.R s9-1-1-transcriptional-and-translational-regulation_data_preprocessing.rmd
Rscript ../run_rmd.R s9-1-2-transcriptional-and-translational-regulation_VHL.rmd
Rscript ../run_rmd.R s9-1-3-transcriptional-and-translational-regulation_EIF4E2.rmd
Rscript ../run_rmd.R s9-1-4-transcriptional-and-translational-regulation_intersection-of-VHL-and-mTOR-v2.rmd

Rscript ../run_rmd.R s9-2-1-alt-TSS-and-HIF.rmd
Rscript ../run_rmd.R s9-2-2-alt-TSS-translation-v4.rmd
Rscript ../run_rmd.R s9-3-1-mTOR-and-TSS-isoforms.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
