s1-2 Preprocessing of fastq files
================
Yoichiro Sugimoto
27 April, 2022

  - [Strategy](#strategy)
  - [Set up](#set-up)
  - [Step 3-2. Copy the first base to read
    name](#step-3-2.-copy-the-first-base-to-read-name)
  - [Step4-0. Collapse technical
    replicates](#step4-0.-collapse-technical-replicates)
  - [Step 4. Map reads to rRNAs and ERCC
    RNAs](#step-4.-map-reads-to-rrnas-and-ercc-rnas)
  - [Session Info](#session-info)

# Strategy

The fastq files will be processed before the alignments to the human
genome. The following steps will be performed:

  - Step 1: Extract UMIs from reads
  - Step 2: Demultiplex reads
  - Step 3: Process index
      - Step 3-1
          - Remove constant regions
          - Demultplex library
      - Step 3-2
          - Copy the first base of reads to read.name
  - Step 4: Map reads to rRNAs and ERCC RNAs
      - Collapse technical replicates
      - ERCC spike-in mapped reads will be used later

This script will perform Step 3-2 and 4.

# Set up

``` r
processors <- 24

library("ShortRead")
library("stringdist") 
library("matrixStats")

temp <- sapply(list.files("../functions", full.names = TRUE), source)
```

``` r
annot.dir <- file.path("../../annotation/")
rRNA.ercc.annot.dir <- file.path(annot.dir, "rRNA_ERCC_annotation")
rRNA.ercc.index.dir <- file.path(rRNA.ercc.annot.dir, "rRNA_ERCC_indices")

results.dir <- file.path("../../results")
processed.fq.dir <- file.path(results.dir, "s1-processed_fastq")
processed.fq.step3.1.dir <- file.path(processed.fq.dir, "s1-1-Step3-1")
processed.fq.step3.2.dir <- file.path(processed.fq.dir, "s1-1-Step3-2")
processed.fq.step4.0.dir <- file.path(processed.fq.dir, "s1-1-Step4-0") 
processed.fq.step4.dir <- file.path(processed.fq.dir, "s1-1-Step4")
processed.fq.step4.2.dir <- file.path(processed.fq.dir, "s1-1-Step4-2-rRNA-ERCC")


create.dirs(
    dirs = c(
        results.dir,
        processed.fq.step3.2.dir,
        processed.fq.step4.0.dir,
        processed.fq.step4.dir,
        processed.fq.step4.2.dir
    )
)

sample.dt <- file.path(
    "../../data/sample_data/processed_sample_file.csv"
) %>% fread

read.process.sample.dt <- file.path(
    "../../data/sample_data/to_process_reads_sample_file.csv"
) %>% fread

## Don't process incomplete data
read.process.sample.dt <- read.process.sample.dt[
    no_technical_rep_sample_name %in% sample.dt[, sample_name]
]
```

# Step 3-2. Copy the first base to read name

This function has a very high memory consumption.

``` r
runAddFirstBaseToFqEntryName.rscript.cmd <- function(analyzed.sample.name, processed.fq.step3.1.dir, processed.fq.step3.2.dir){

    addFirstBaseToFqEntryName.rscript <- file.path("./functions/addFirstBaseToFqEntryName.R")

    run.addFirstBaseToFqEntryName.rscript.cmd <- paste(
        "Rscript",
        addFirstBaseToFqEntryName.rscript,
        "-s", analyzed.sample.name,
        "-i", processed.fq.step3.1.dir,
        "-o", processed.fq.step3.2.dir
    )

    system.cat(run.addFirstBaseToFqEntryName.rscript.cmd)
    
    return()
}


temp <- mclapply(
    read.process.sample.dt[, sample_name_with_lane],
    runAddFirstBaseToFqEntryName.rscript.cmd,
    processed.fq.step3.1.dir = processed.fq.step3.1.dir,
    processed.fq.step3.2.dir = processed.fq.step3.2.dir,
    mc.cores = round(processors / 4)
)
```

# Step4-0. Collapse technical replicates

``` r
collapseTechnicalReplicates <- function(no.techr.sample.name, read.process.sample.dt, processed.fq.step3.2.dir, processed.fq.step4.0.dir){
    input.files <- file.path(
        processed.fq.step3.2.dir,
        paste0(read.process.sample.dt[
            no_technical_rep_sample_name == no.techr.sample.name,
            sample_name_with_lane
        ], "_R1.fastq.gz")
    )

    output.file <- file.path(
        processed.fq.step4.0.dir,
        paste0(no.techr.sample.name, "_R1.fastq.gz")
    )

    cat.command <- paste(
        c("cat", input.files, ">", output.file),
        collapse = " "
    ) %T>%
        system.cat

    gsub("_R1.fastq.gz", "_R2.fastq.gz", cat.command) %>%
        system.cat
    
    return()
}

## Sanity check
if(all(
    sort(sample.dt[, sample_name]) ==
    sort(read.process.sample.dt[
        !duplicated(no_technical_rep_sample_name), no_technical_rep_sample_name
    ])
)){"OK"} else {stop()}
```

    ## [1] "OK"

``` r
temp <- mclapply(
    sample.dt[, sample_name],
    collapseTechnicalReplicates,
    read.process.sample.dt = read.process.sample.dt,
    processed.fq.step3.2.dir = processed.fq.step3.2.dir,
    processed.fq.step4.0.dir = processed.fq.step4.0.dir,
    mc.cores = round(processors / 2)
)
```

# Step 4. Map reads to rRNAs and ERCC RNAs

``` r
removeRibosomalRNA <- function(analyzed.sample.name, rRNA.ercc.index.dir, processed.fq.step4.0.dir, processed.fq.step4.dir, processed.fq.step4.2.dir, processors){

    step4.in.files <- file.path(
        processed.fq.step4.0.dir,
        paste0(analyzed.sample.name, c("_R1.fastq.gz", "_R2.fastq.gz"))
    )
    
    step4.out.files <- file.path(
        processed.fq.step4.dir,
        gsub("_R1", "", basename(step4.in.files[1]))
    )

    step4.cmd <- paste(
        "bowtie2",
        "-p", processors,
        "--un-conc-gz", step4.out.files,
        "-N", 1,
        "-x", paste0(rRNA.ercc.index.dir, "/rRNA"),
        "-1", step4.in.files[1],
        "-2", step4.in.files[2],
        "| samtools view -@", processors, "-Su /dev/stdin",
        "| samtools sort -@", processors, "-T", analyzed.sample.name, "-",
        ">", file.path(processed.fq.step4.2.dir, paste0(analyzed.sample.name, ".sorted.bam"))
    )

    ## cat(analyzed.sample.name, sep = "\n\n")
    step4.out <- system.cat(step4.cmd)
    return()
}

temp <- lapply(
    sample.dt[, sample_name],
    removeRibosomalRNA,
    rRNA.ercc.index.dir = rRNA.ercc.index.dir,
    processed.fq.step4.0.dir = processed.fq.step4.0.dir,
    processed.fq.step4.dir = processed.fq.step4.dir,
    processed.fq.step4.2.dir = processed.fq.step4.2.dir,
    processors = processors
)
```

# Session Info

``` r
sessionInfo()
```

    ## R version 4.0.0 (2020-04-24)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/five_prime_seq_for_VHL_loss_v0.1.1/lib/libopenblasp-r0.3.9.so
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
    ##  [1] knitr_1.28                  stringr_1.4.0              
    ##  [3] magrittr_1.5                data.table_1.12.8          
    ##  [5] dplyr_1.0.0                 khroma_1.3.0               
    ##  [7] ggplot2_3.3.1               stringdist_0.9.5.5         
    ##  [9] ShortRead_1.46.0            GenomicAlignments_1.24.0   
    ## [11] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
    ## [13] matrixStats_0.56.0          Biobase_2.48.0             
    ## [15] Rsamtools_2.4.0             GenomicRanges_1.40.0       
    ## [17] GenomeInfoDb_1.24.0         Biostrings_2.56.0          
    ## [19] XVector_0.28.0              IRanges_2.22.1             
    ## [21] S4Vectors_0.26.0            BiocParallel_1.22.0        
    ## [23] BiocGenerics_0.34.0         rmarkdown_2.2              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0       xfun_0.14              purrr_0.3.4           
    ##  [4] lattice_0.20-41        colorspace_1.4-1       vctrs_0.3.1           
    ##  [7] generics_0.0.2         htmltools_0.4.0        yaml_2.2.1            
    ## [10] rlang_0.4.6            pillar_1.4.4           withr_2.2.0           
    ## [13] glue_1.4.1             RColorBrewer_1.1-2     jpeg_0.1-8.1          
    ## [16] GenomeInfoDbData_1.2.3 lifecycle_0.2.0        zlibbioc_1.34.0       
    ## [19] munsell_0.5.0          gtable_0.3.0           hwriter_1.3.2         
    ## [22] evaluate_0.14          latticeExtra_0.6-29    Rcpp_1.0.4.6          
    ## [25] scales_1.1.1           png_0.1-7              digest_0.6.25         
    ## [28] stringi_1.4.6          grid_4.0.0             tools_4.0.0           
    ## [31] bitops_1.0-6           RCurl_1.98-1.2         tibble_3.0.1          
    ## [34] crayon_1.3.4           pkgconfig_2.0.3        ellipsis_0.3.1        
    ## [37] Matrix_1.2-18          R6_2.4.1               compiler_4.0.0
