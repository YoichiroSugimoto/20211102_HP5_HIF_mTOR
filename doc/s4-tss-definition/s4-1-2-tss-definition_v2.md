---
title: "s4-1-2 Definition of tss (count per cluster)"
author: "Yoichiro Sugimoto"
date: "13 November, 2021"
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
		
The count per TSS will be calculated.



```r
knitr::opts_chunk$set(collapse = TRUE)
```


```r

## Data analysis
library("data.table")
library("magrittr")
library("stringr")
## Bioconductor
library("rtracklayer")
library("GenomicFeatures")
library("GenomicAlignments")
library("Rsubread")
## Data visualization
library("ggplot2")
## Parallelization
library("parallel")

analysis.variable.dt <- fread("../../data/sample_data/analysis_variables.csv")

## chrom.size.file <- "../../annotation/hg38_annotation/star_indices/chrNameLength.txt"
## chrom.size.dt <- fread(chrom.size.file)

## Specify the number of CPUs to be used
processors <- 8

temp <- sapply(list.files("../functions", full.names = TRUE), source)

sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

results.dir <- file.path("../../results")

## Input data
s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
## star.aligned.bam.dir <- file.path(s2.alignment.dir, "s2-1-b-star-aligned_bam")
s2.2.processed.bam.dir <-  file.path(s2.alignment.dir, "s2-2-processed-data")

tss.bam.dir <- ifelse(
    analysis.variable.dt[variable_name == "umi_dedup", selected_variable] == "dedup",
    file.path(s2.2.processed.bam.dir, "s2-2-2-dedup-tss-bam"),
    file.path(s2.2.processed.bam.dir, "s2-2-1-tss-bam")
)

## Define output directory
s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.4.dpi.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-4-tss-by-dpi")
s4.1.4.1.dpi.in.bed.dir <- file.path(s4.1.4.dpi.tss.dir, "dpi-input-bed")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")
s4.1.7.count.per.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-7-count-per-tss") 

create.dirs(c(
    s4.1.7.count.per.tss.dir
))
```


# Calculate count per TSS of all sequencing data



```r

sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]

filtered.cluster.with.quantile.file <- file.path(
    s4.1.6.filtered.tss.dir,
    "filtered-tss-with-quantile.csv"
)

filtered.tss.with.quantile.dt <- fread(filtered.cluster.with.quantile.file)

count.per.tss.dt <- countPerRangeFromBams(
    list.files(tss.bam.dir, pattern = "tss.bam$", full.names = TRUE),
    gr = makeGRangesFromDataFrame(filtered.tss.with.quantile.dt, keep.extra.columns = TRUE),
    out.file.name = file.path(s4.1.7.count.per.tss.dir, "count-per-confident-tss.csv"),
    sample.names = sample.names,
    processors = processors,
    range.name = "tss_name",
    bam.file.prefix = ifelse(
        analysis.variable.dt[variable_name == "umi_dedup", selected_variable] == "dedup",
        ".tss.dedup.tss.bam$",
        ".tss.bam$"
    )
)

```

# Calculate count per TSS of FANTOM5 CAGE data



```r


countPerRangeFromBeds <- function(input.bed.files, gr, out.file.name, processors, range.name = "tss_name", bed.file.sample.indicator = "CNhs[[:digit:]]*", sample.names = NULL){
    ## This function count the number of TSS bam mapped within ranges specified by genomic ranges
    ## input.bam.files: vector of input bam files
    ## ge: GRanges of range
    ## out.file.name: output count table filename
    ## processors: the number of processors to be used
    ## range.name: the names of range
    ## bed.file.postfix: input bedfile is expected to have the name of sample.name + bed.file.postfix
    ## sample.names: if sample.names are specified here, the order of column will be sorted in this order
    require("GenomicAlignments")
    
    createCountTable <- function(input.bed.file, gr, range.name, bed.file.sample.indicator){
        bed.dt <- fread(cmd = paste0("zcat ", input.bed.file))
        setnames(
            bed.dt,
            old = paste0("V", 1:6),
            new = c("chr", "start", "end", "cluster_name", "score", "strand")
        )
        bed.dt <- bed.dt[rep(1:nrow(bed.dt), times = score)]
        
        ol.count <- countOverlaps(gr, makeGRangesFromDataFrame(bed.dt))
        ol.count.dt <- data.table(V1 = mcols(gr)[, range.name], count = ol.count)

        ol.count.dt <- ol.count.dt[, .(count = sum(count)), by = V1]
        setnames(
            ol.count.dt,
            old = c("V1", "count"),
            new = c(range.name, str_extract(basename(input.bed.file), bed.file.sample.indicator))
        )
        setkeyv(ol.count.dt, range.name)
        return(ol.count.dt)
    }

    input.bed.files <- input.bed.files[
        grep(bed.file.sample.indicator, input.bed.files)
    ]

    if(all(grep(bed.file.sample.indicator, input.bed.files))){
        "OK"
    } else {
        print("The following files does not have the file name indicator")
        print(input.bed.files[!grepl(bed.file.sample.indicator, input.bed.files)])
    }
    
    count.dts <- mclapply(
        input.bed.files,
        createCountTable,
        gr = gr,
        range.name = range.name,
        bed.file.sample.indicator = bed.file.sample.indicator,
        mc.cores = processors
    )

    count.dt <- Reduce(function(...) merge(..., by = range.name), count.dts)
    setkeyv(count.dt, range.name)

    if(!is.null(sample.names)){ 
        count.dt <- count.dt[, c(range.name, sample.names), with = FALSE]
    } else {"No columnname sorting"}
    
    fwrite(count.dt, file = out.file.name)
    
    return(count.dt)
}

fantom.cage.count.per.tss.dt <- list.files(
    s4.1.4.1.dpi.in.bed.dir,
    pattern = "ctss.bed.gz$",
    full.names = TRUE
) %>%
    {.[!grepl("^total_", basename(.))]} %>%
    {countPerRangeFromBeds(
         input.bed.files = .,
         gr = makeGRangesFromDataFrame(filtered.tss.with.quantile.dt, keep.extra.columns = TRUE),
         out.file.name = file.path(s4.1.7.count.per.tss.dir, "FANTOM-CAGE-count-per-confident-tss.csv"),
         processors = processors,
         range.name = "tss_name",
         bed.file.sample.indicator = "CNhs[[:digit:]]*"
     )}

```


 
# Session information



```r

sessionInfo()
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/run_dpi_v2/lib/libopenblasp-r0.3.10.so
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
##  [1] knitr_1.29                  dplyr_1.0.0                
##  [3] khroma_1.4.0                ggplot2_3.3.2              
##  [5] Rsubread_2.0.0              GenomicAlignments_1.22.0   
##  [7] Rsamtools_2.2.0             Biostrings_2.54.0          
##  [9] XVector_0.26.0              SummarizedExperiment_1.16.0
## [11] DelayedArray_0.12.0         BiocParallel_1.20.0        
## [13] matrixStats_0.56.0          GenomicFeatures_1.38.0     
## [15] AnnotationDbi_1.48.0        Biobase_2.46.0             
## [17] rtracklayer_1.46.0          GenomicRanges_1.38.0       
## [19] GenomeInfoDb_1.22.0         IRanges_2.20.0             
## [21] S4Vectors_0.24.0            BiocGenerics_0.32.0        
## [23] stringr_1.4.0               magrittr_1.5               
## [25] data.table_1.12.8           rmarkdown_2.3              
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5             lattice_0.20-41        prettyunits_1.1.1     
##  [4] assertthat_0.2.1       digest_0.6.25          BiocFileCache_1.10.0  
##  [7] R6_2.4.1               RSQLite_2.2.0          evaluate_0.14         
## [10] httr_1.4.1             pillar_1.4.6           zlibbioc_1.32.0       
## [13] rlang_0.4.7            progress_1.2.2         curl_4.3              
## [16] blob_1.2.1             Matrix_1.2-18          munsell_0.5.0         
## [19] RCurl_1.98-1.2         bit_1.1-15.2           biomaRt_2.42.0        
## [22] compiler_3.6.3         xfun_0.15              pkgconfig_2.0.3       
## [25] askpass_1.1            htmltools_0.5.0        openssl_1.4.2         
## [28] tidyselect_1.1.0       tibble_3.0.3           GenomeInfoDbData_1.2.2
## [31] XML_3.99-0.3           withr_2.2.0            crayon_1.3.4          
## [34] dbplyr_1.4.4           bitops_1.0-6           rappdirs_0.3.1        
## [37] grid_3.6.3             gtable_0.3.0           lifecycle_0.2.0       
## [40] DBI_1.1.0              scales_1.1.1           stringi_1.4.6         
## [43] ellipsis_0.3.1         vctrs_0.3.2            generics_0.0.2        
## [46] tools_3.6.3            bit64_0.9-7.1          glue_1.4.1            
## [49] purrr_0.3.4            hms_0.5.3              yaml_2.2.1            
## [52] colorspace_1.4-1       memoise_1.1.0
```
