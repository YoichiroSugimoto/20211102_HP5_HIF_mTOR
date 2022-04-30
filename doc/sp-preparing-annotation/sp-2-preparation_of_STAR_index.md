---
title: "sp-2 Preparation of STAR index"
author: "Yoichiro Sugimoto"
date: "27 April, 2022"
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

The sequences of genome are downloaded, and exported as
a `fasta` format.
The retrieved sequences are then indexed to map high-throughput DNA sequencing
reads. 
To do this, [STAR](https://github.com/alexdobin/STAR) will be used.

# Setup



```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

### Capture `system` outputs
system <- function(...) {
    stopifnot(!any(names(list(...)) %in% "intern"))
    result <- base::system(..., intern = TRUE)
    cat(paste0(result, "\n"))
}
```



```r
file.prefix <- "hg38"
### human genome by BSgenome
bs.genome.ver <- "BSgenome.Hsapiens.UCSC.hg38"
library(bs.genome.ver, character.only = TRUE)

processors <- 8
```

## The input and output files and the directories



```r
## .libPaths("/well/ratcliff/data/yoichiro/software/R_packages/library")
annot.dir <- file.path("../../annotation/")
ref.seq.dir <- file.path(annot.dir, "hg38_annotation/ref_sequences/")
star.index.dir <- file.path(annot.dir, "hg38_annotation/star_indices/") 
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")

create.dir <- function(dir.name){
    if(dir.exists(dir.name) == FALSE) {
        dir.create(dir.name)
    }
}
create.dir(ref.seq.dir)
create.dir(star.index.dir)
```



## Load pre-processed data

The annotation data pre-processed by "i Parse gene annotation" script
is loaded here.



```r
all.tx.gtf <- file.path(
    annot.ps.dir,
    "all-transcript.gtf"
)
```

# Genome

Follow STAR manual

"It is strongly recommended to include major chromosomes (e.g., for human chr1-22,chrX,chrY,chrM,) as well as un-placed and un-localized scaffolds. Typically, un-placed/un-localized scaffolds add just a few MegaBases to the genome length, however, a substantial number of reads may map to ribosomal RNA (rRNA) repeats on these scaffolds. These reads would be reported as unmapped if the scaffolds are not included in the genome, or, even worse, may be aligned to wrong loci on the chromosomes. Generally, patches and alternative haplotypes should not be included in the genome."

In BSgenome file, the following prefixes were found (explanation was taken from [NCBI](https://www.ncbi.nlm.nih.gov/grc/help/patches/) and [UCSC](http://hgdownload.soe.ucsc.edu/gbdb/hg38/html/description.html)):

 - alt
	 - "Haplotype chromosome, unplaced contig and unlocalized contig names now include their NCBI accession number (e.g., chr6_GL000256v2_alt)"
	 - "Haplotype chromosome names consist of the chromosome number, followed by the NCBI accession number, followed by alt"
 - fix
	 - "FIX patches: Fix patches represent changes to existing assembly sequences."
 - chrUn
	 - "Unplaced contig names (contigs whose associated chromosome is not known) consist of "chrUn" followed by the NCBI accession number"
 - random
	 - "Unlocalized contig names consist of the chromosome number, followed by the NCBI accession number, followed by random"

So, I exclude chromosomes where the names containing the postfix of _alt and _fix.


```r
genome.fa.file <- file.path(ref.seq.dir, paste0(file.prefix, "_genome.fa"))

## The original function extracted standard chromosomes only but here I include all chrosomesomes
## eval(substitute(
##     genome.seq <- sapply(chromosomes, function(x){bsgenome.ver[[x]]}),
##     list(bsgenome.ver = as.name(bs.genome.ver))
## ))

hg38.all.chrs <- names(getBSgenome(BSgenome.Hsapiens.UCSC.hg38))

hg38.chrs <- hg38.all.chrs %>%
    .[!grepl("_alt$", .)] %>%
    .[!grepl("_fix$", .)]

eval(substitute(
    genome.seq <- sapply(hg38.chrs, function(x){bsgenome.ver[[x]]}),
    list(bsgenome.ver = as.name(bs.genome.ver))
))

export(genome.seq, con = genome.fa.file, format = "fasta")
```

# Index retrieved sequences

Here, index for `star` is generated. 


```r
## Build STAR index
buildSTARIndex <- function(star.cmd = "STAR", star.index.dir, genome.fa.file, all.tx.gtf){
    cmd <- paste(star.cmd,
                 "--runThreadN", processors,
                 "--runMode genomeGenerate",
                 "--genomeDir", star.index.dir,
                 "--genomeFastaFiles", genome.fa.file,
                 "--sjdbGTFfile", all.tx.gtf,
                 "--sjdbOverhang 99")
	system(cmd) 
}

buildSTARIndex(star.cmd = "STAR", star.index.dir, genome.fa.file, all.tx.gtf)
```

```
## Apr 27 14:46:42 ..... started STAR run
##  Apr 27 14:46:42 ... starting to generate Genome files
##  Apr 27 14:47:44 ... starting to sort Suffix Array. This may take a long time...
##  Apr 27 14:48:00 ... sorting Suffix Array chunks and saving them to disk...
##  Apr 27 15:38:03 ... loading chunks from disk, packing SA...
##  Apr 27 15:39:03 ... finished generating suffix array
##  Apr 27 15:39:03 ... generating Suffix Array index
##  Apr 27 15:42:17 ... completed Suffix Array index
##  Apr 27 15:42:17 ..... processing annotations GTF
##  Apr 27 15:42:31 ..... inserting junctions into the genome indices
##  Apr 27 15:46:02 ... writing Genome to disk ...
##  Apr 27 15:46:06 ... writing Suffix Array to disk ...
##  Apr 27 15:46:20 ... writing SAindex to disk
##  Apr 27 15:46:21 ..... finished successfully
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.3 BSgenome_1.56.0                  
##  [3] rtracklayer_1.48.0                Biostrings_2.56.0                
##  [5] XVector_0.28.0                    GenomicRanges_1.40.0             
##  [7] GenomeInfoDb_1.24.0               IRanges_2.22.1                   
##  [9] S4Vectors_0.26.0                  BiocGenerics_0.34.0              
## [11] knitr_1.28                        stringr_1.4.0                    
## [13] magrittr_1.5                      data.table_1.12.8                
## [15] dplyr_1.0.0                       khroma_1.3.0                     
## [17] ggplot2_3.3.1                     rmarkdown_2.2                    
## 
## loaded via a namespace (and not attached):
##  [1] SummarizedExperiment_1.18.1 tidyselect_1.1.0           
##  [3] xfun_0.14                   purrr_0.3.4                
##  [5] lattice_0.20-41             colorspace_1.4-1           
##  [7] vctrs_0.3.1                 generics_0.0.2             
##  [9] htmltools_0.4.0             yaml_2.2.1                 
## [11] XML_3.99-0.3                rlang_0.4.10               
## [13] pillar_1.4.4                glue_1.4.1                 
## [15] withr_2.4.1                 BiocParallel_1.22.0        
## [17] matrixStats_0.56.0          GenomeInfoDbData_1.2.3     
## [19] lifecycle_0.2.0             zlibbioc_1.34.0            
## [21] munsell_0.5.0               gtable_0.3.0               
## [23] evaluate_0.14               Biobase_2.48.0             
## [25] Rcpp_1.0.4.6                scales_1.1.1               
## [27] DelayedArray_0.14.0         Rsamtools_2.4.0            
## [29] digest_0.6.25               stringi_1.4.6              
## [31] grid_4.0.0                  tools_4.0.0                
## [33] bitops_1.0-6                RCurl_1.98-1.2             
## [35] tibble_3.0.1                crayon_1.3.4               
## [37] pkgconfig_2.0.3             ellipsis_0.3.1             
## [39] Matrix_1.2-18               R6_2.4.1                   
## [41] GenomicAlignments_1.24.0    compiler_4.0.0
```
