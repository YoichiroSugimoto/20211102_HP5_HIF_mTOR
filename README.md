# Introduction

This repository contains the pipeline for **"Isoform-resolved mRNA profiling of ribosome load defines interplay of HIF and mTOR dysregulation in kidney cancer"** (Sugimoto and Ratcliffe).

This pipeline can generate most of plots and results presented in the article using the following sequence data deposited in ArrayExpress as input:

- HP5: E-MTAB-10689
- 5′ end-Seq of total mRNAs: E-MTAB-10688


# Usage

The script to run the entire pipeline is located at `master_scripts` directory. The output of the analyses is stored in `doc` directory.

# Prerequisite

Two `conda` environments summarised in `conda_environment` must be created before the analysis. A few additional packages may be required for some of the analyses.

## Required data

- Raw sequence data (available from ArrayExpress, HP5: E-MTAB-10689, 5′ end-Seq of total mRNAs: E-MTAB-10688)
- ERCC RNA sequence: the sequences downloaded from NIST (https://www-s.nist.gov/srmors/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v1.fasta) should be stored at `data/ERCC_sequence`. The sequence file is also available from this github repository.
- Human CAGE data from FANTOM5

## Required software
In addition to the software included in the `conda` environments, the following software is required to run this pipeline.

- [paraclu](http://cbrc3.cbrc.jp/~martin/paraclu/)
- [dpi](https://github.com/hkawaji/dpi1)
    - I slightly rewrote some of the script of dpi so that the software can be run on `Slurm`. The modified version is available from the following link (https://github.com/YoichiroSugimoto/dpi).

# Analysis overview

## Structure

The analyses consist of the following steps:

- sp-0 Annotation preparation
- s1 sequence read preprocessing
- s2 Sequence read alignment
- s3 Alignment statistics
- s4 TSS definition
- s5 TSS validation
- s6 Analysis of differential TSS usage by VHL loss or hypoxia
- s7 HIF binding sites definition
- s8 Analysis of translation
- s9 Integrative analysis
- s10 Additional analysis

## **sp-0 Annotation preparation**

Annotation used throughout the analyses is prepared.


## **s1 Sequence read preprocessing**

`fastq` files of 5' end-Seq are preprocessed for the downstream analysis. This includes the adapter trimming, UMI extraction, and removal of reads from cytoplasmic and mitochondrial rRNAs.

(Note) This script is very slow taking ~24h with 18 cores and 140G memory.
After running this script, make sure that the memories had not been exhausted during the process by checking the `stdout`.


## **s2 Sequence read alignment**

The preprocessed sequence reads are mapped to the human genome.

(Note) This script tries to remove `reference genome` from the memory before running `star`. This is a precautious measure not to use any unexpected genome from the memory and I found it was necessary at some setting. However, this generates a harmless error `EXITING: Did not find the genome in memory, did not remove any genomes from shared memory` if there is no genome in the memory and this message will be included in the output.


## **s3 Alignment statistics**


## **s4 TSS definition**

TSS are defined based on 5' end-Seq data. In addition, a transcript is assigned to a TSS, and various features of the transcripts are analysed.


## **s5 TSS validation**

The 5' end-Seq defined TSS are validated by comparing annotated TSS by RefSeq and GENCODE.


## **s6 Analysis of differential TSS usage by VHL loss or hypoxia**

The differential TSS usage upon VHL loss or hypoxia are examined.


## **s7 HIF binding sites definition**

HIF binding sites and HIF2A/HIF1A binding ratio are examined.


## **s8 Analysis of translation**

Changes in mean ribosome load (i.e. absolute translational efficiency) are calculated.

In addition, validation of the HP5 method is performed.

Finally, the effect of mTOR inhibition on translational efficiency is evaluated using HP5.

## **s9 Integrative analysis**

Most of analyses used for the publication are performed in this section.

## **s10 Additional analysis**
