---
title: "sp-3 Preparation of bowtie2 index for the alignments to rRNA and ERCC"
author: "Yoichiro Sugimoto"
date: "28 February, 2022"
vignette: >
  %\VignetteIndexEntry{Bioconductor style for PDF documents}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
   html_document:
     highlight: haddock
     toc: yes
     toc_depth: 2
     keep_md: yes
     fig_width: 5
     fig_height: 5
---

# Overview

The rRNA and [ERCC spike-in RNA](https://www-s.nist.gov/srmors/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v1.fasta) sequences are downloaded.

The retrieved sequence is then indexed to map sequence reads. 
To do this, bowtie2 will be used.

# Setup

## The input and output files and the directories


```r
library("rentrez")
library("rtracklayer")
library("Biostrings")

temp <- sapply(list.files("../functions", full.names = TRUE), source)

processors <- 8
```


```r
annot.dir <- file.path("../../annotation/")
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")

rRNA.ercc.annot.dir <- file.path(annot.dir, "rRNA_ERCC_annotation")
rRNA.ercc.index.dir <- file.path(rRNA.ercc.annot.dir, "rRNA_ERCC_indices")

create.dirs(c(
    rRNA.ercc.annot.dir,
    rRNA.ercc.index.dir
))
```


# Download human rRNA sequences and generate rRNA and ERCC bowtie2 index



```r
## Normal rRNAs
rRNA.ids <- c("NR_023363.1","NR_046235.1") # NCBI rRNA 5S, pre-rRNA 45S: NR_046235.1
## Since 5.8S, 28S, 18S is included in 45S, I did not include "NR_003285.2","NR_003287.2","NR_003286.2"
rRNA.fa <- entrez_fetch(db="nucleotide", id=rRNA.ids, rettype="fasta")
rRNA.fa.file <- file.path(
    rRNA.ercc.annot.dir,
    paste0(
        "rRNA_",
        format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
        ".fa"
    ))
write(rRNA.fa, file = rRNA.fa.file)

## mitochondria rRNAs
primary.tx.gtf <- list.files(
    file.path(annot.dir, "hg38_annotation/processed_data"),
    pattern = "primary_transcript_[[:digit:]]*_[[:digit:]]*_[[:digit:]]*\\.gtf$",
    full.names = TRUE
) 

primary.tx.annotation.gr <-import(primary.tx.gtf)

bs.genome.ver <- "BSgenome.Hsapiens.UCSC.hg38"
library(bs.genome.ver, character.only = TRUE)
```

```
## Loading required package: BSgenome
```

```r
mt.rRNA.ids <- c("ENSG00000211459", "ENSG00000210082") #MT-RNR1, MT-RNR2
mt.rRNA.gr <- primary.tx.annotation.gr[primary.tx.annotation.gr$gene_id %in% mt.rRNA.ids]

mt.rRNA.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, mt.rRNA.gr)
names(mt.rRNA.seq) <- mt.rRNA.ids

export(mt.rRNA.seq, con = rRNA.fa.file, format = "fasta", append = TRUE)

ercc.fa.file <- file.path(
    "../../data/ERCC_sequence/SRM2374_putative_T7_products_NoPolyA_v1.fasta.txt"
)

ercc.seq <- readDNAStringSet(ercc.fa.file, format = "fasta")

export(ercc.seq, con = rRNA.fa.file, format = "fasta", append = TRUE)

## bowtie2 indexing

index.cmd <- paste(
    "bowtie2-build",
    rRNA.fa.file,
    file.path(rRNA.ercc.index.dir, "rRNA")
)

index.out <- system.cat(index.cmd)
```


# Session information


```r
sessionInfo()
```

```
## R version 4.0.0 (2020-04-24)
## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/five_prime_seq_for_VHL_loss_v0.2.1/lib/libopenblasp-r0.3.10.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.3 BSgenome_1.56.0                  
##  [3] knitr_1.28                        stringr_1.4.0                    
##  [5] magrittr_1.5                      data.table_1.12.8                
##  [7] dplyr_1.0.0                       khroma_1.3.0                     
##  [9] ggplot2_3.3.1                     Biostrings_2.56.0                
## [11] XVector_0.28.0                    rtracklayer_1.48.0               
## [13] GenomicRanges_1.40.0              GenomeInfoDb_1.24.0              
## [15] IRanges_2.22.1                    S4Vectors_0.26.0                 
## [17] BiocGenerics_0.34.0               rentrez_1.2.2                    
## [19] rmarkdown_2.2                    
## 
## loaded via a namespace (and not attached):
##  [1] SummarizedExperiment_1.18.1 tidyselect_1.1.0           
##  [3] xfun_0.14                   purrr_0.3.4                
##  [5] lattice_0.20-41             colorspace_1.4-1           
##  [7] vctrs_0.3.1                 generics_0.0.2             
##  [9] htmltools_0.4.0             yaml_2.2.1                 
## [11] XML_3.99-0.3                rlang_0.4.10               
## [13] pillar_1.4.4                withr_2.4.1                
## [15] glue_1.4.1                  BiocParallel_1.22.0        
## [17] matrixStats_0.56.0          GenomeInfoDbData_1.2.3     
## [19] lifecycle_0.2.0             zlibbioc_1.34.0            
## [21] munsell_0.5.0               gtable_0.3.0               
## [23] evaluate_0.14               Biobase_2.48.0             
## [25] curl_4.3                    Rcpp_1.0.4.6               
## [27] scales_1.1.1                DelayedArray_0.14.0        
## [29] jsonlite_1.7.2              Rsamtools_2.4.0            
## [31] digest_0.6.25               stringi_1.4.6              
## [33] grid_4.0.0                  tools_4.0.0                
## [35] bitops_1.0-6                RCurl_1.98-1.2             
## [37] tibble_3.0.1                crayon_1.3.4               
## [39] pkgconfig_2.0.3             ellipsis_0.3.1             
## [41] Matrix_1.2-18               httr_1.4.2                 
## [43] R6_2.4.1                    GenomicAlignments_1.24.0   
## [45] compiler_4.0.0
```
