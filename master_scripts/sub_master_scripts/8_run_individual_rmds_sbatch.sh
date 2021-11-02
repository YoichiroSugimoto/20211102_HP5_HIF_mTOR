#!/bin/bash
#SBATCH --job-name=8_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --time=3-00:00:0
#SBATCH --mem=192G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch --dependency=afterok:25477895 /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/master_scripts/sub_master_scripts/8_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210219_HP5_VHL_mTOR/R/s8-analysis-of-translation

Rscript ../run_rmd.R s8-1-1-gene-level-translation-analysis.rmd
Rscript ../run_rmd.R s8-1-2-tx-level-translation-analysis.rmd
Rscript ../run_rmd.R s8-2-identification-of-dte-isoforms.rmd
Rscript ../run_rmd.R s8-3-0-validation-of-method_data-preprocessing.rmd
Rscript ../run_rmd.R s8-3-1-1-validation-of-method_single-condition-MRL-prediction-v3.rmd
Rscript ../run_rmd.R s8-3-1-2-validation-of-method_single-condition-isoform-comparison.rmd
Rscript ../run_rmd.R s8-3-2-1-analysis-of-mTOR-dependent-translational-regulation-1.rmd
Rscript ../run_rmd.R s8-3-2-2-analysis-of-mTOR-dependent-translational-regulation-2.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
