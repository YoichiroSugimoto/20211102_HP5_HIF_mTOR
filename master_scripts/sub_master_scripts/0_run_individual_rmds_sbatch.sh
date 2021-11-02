#!/bin/bash
#SBATCH --job-name=0_five_prime_seq_for_VHL_loss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/master_scripts/sub_master_scripts/0_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate five_prime_seq_for_VHL_loss_v0.2.1

## Preprocessing annotations
cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/R/sp-preparing-annotation

Rscript ../run_rmd.R sp-1-parse_gene_annotation.rmd
Rscript ../run_rmd.R sp-2-preparation_of_STAR_index.rmd
Rscript ../run_rmd.R sp-3-preparation_of_bowtie2_index.rmd
Rscript ../run_rmd.R sp-4-preparation_of_KEGG_mapping.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
