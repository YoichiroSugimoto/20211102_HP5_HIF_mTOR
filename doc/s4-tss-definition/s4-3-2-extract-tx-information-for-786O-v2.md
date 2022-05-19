s4-3-2 Extract transcript information for 786O
================
Yoichiro Sugimoto
18 May, 2022

  - [Overview](#overview)
  - [Define variables](#define-variables)
  - [Calculate read count for each position of
    TSS](#calculate-read-count-for-each-position-of-tss)
  - [Extract transcript meta
    information](#extract-transcript-meta-information)
  - [Session information](#session-information)

# Overview

``` r
library("rtracklayer")
library("GenomicFeatures")
library("GenomicAlignments")

## Specify the number of CPUs to be used
processors <- 16

temp <- sapply(list.files("../functions", full.names = TRUE), source)
source("./functions/extractTxMetaInfo.R")

sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

annot.dir <- normalizePath(file.path("../../annotation/"))
ref.seq.dir <- file.path(annot.dir, "hg38_annotation/ref_sequences/")
genome.fa.file <- list.files(ref.seq.dir, pattern = "_genome.fa$", full.names = TRUE)

annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)

results.dir <- file.path("../../results")
s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
s2.2.processed.bam.dir <-  file.path(s2.alignment.dir, "s2-2-processed-data")

s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")

s4.2.tx.assignment.dir <- file.path(s4.tss.dir, "s4-2-transcript-assignment")
s4.2.2.tss.tx.map.786O.dir <- file.path(s4.2.tx.assignment.dir, "s4-2-2-tss-transcript-mapping-786O")

s4.3.tx.info.dir <- file.path(s4.tss.dir, "s4-3-transcript-info")
s4.3.2.tx.info.786o.dir <- file.path(s4.3.tx.info.dir, "s4-3-2-transcript-info-for-786O")

create.dirs(c(
    s4.3.tx.info.dir,
    s4.3.2.tx.info.786o.dir
))

sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]
```

# Define variables

``` r
tx.gtf.file <- file.path(
    s4.2.2.tss.tx.map.786O.dir,
    "transcripts-per-TSS-for-786O.gtf"
)

ref.colname <- "tss_name"

tx.info.dir <- s4.3.2.tx.info.786o.dir

tx.meta.data.file <- file.path(
  tx.info.dir,
  "transcript-meta-information-786O-VHL.csv"
)
```

# Calculate read count for each position of TSS

``` r
filtered.tss.with.quantile.dt <- file.path(
  s4.1.6.filtered.tss.dir,
  "filtered-tss-with-quantile.csv"
) %>% fread

analysis.variable.dt <- fread("../../data/sample_data/analysis_variables.csv")
vhl.status <- "VHL"

tss.bam.dir <- ifelse(
    analysis.variable.dt[variable_name == "umi_dedup", selected_variable] == "dedup",
    file.path(s2.2.processed.bam.dir, "s2-2-2-dedup-tss-bam"),
    file.path(s2.2.processed.bam.dir, "s2-2-1-tss-bam")
)

input.bam.files <- list.files(
    tss.bam.dir,
    full.names = TRUE,
    pattern = ifelse(
        analysis.variable.dt[variable_name == "umi_dedup", selected_variable] == "dedup",
        ".tss.dedup.tss.bam$",
        ".tss.bam$"
    )
) %>%
    {grep(paste0("total_786O_", vhl.status, "_HIF1B_N"), ., value = TRUE)}


tss.count.per.pos.dt <- tssWeightByPosition(
    filtered.tss.with.quantile.dt = filtered.tss.with.quantile.dt,
    input.bam.files = input.bam.files
)
```

    ## [1] "The following files will be used:"
    ## [1] "total_786O_VHL_HIF1B_N_1.tss.bam" "total_786O_VHL_HIF1B_N_2.tss.bam"
    ## [3] "total_786O_VHL_HIF1B_N_3.tss.bam" "total_786O_VHL_HIF1B_N_4.tss.bam"

# Extract transcript meta information

``` r
temp <- extractTxMetaInfo(
    tx.gtf.file = tx.gtf.file,
    tss.count.per.pos.dt = tss.count.per.pos.dt,
    ref.colname = ref.colname,
    tx.info.dir = tx.info.dir,
    tx.meta.data.file = tx.meta.data.file,
    processors = processors
)
```

    ## [1] "Export transcript sequences"

    ## Loading required package: BSgenome

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ...

    ## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): The "phase" metadata column contains non-NA values for features of type
    ##   stop_codon. This information was ignored.

    ## Warning in .find_exon_cds(exons, cds): The following transcripts have exons that contain more than one CDS
    ##   (only the first CDS was kept for each exon): ENSG00000104904_1,
    ##   ENSG00000143450_1, ENSG00000180304_1

    ## OK
    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    ## [1] "Analyzing uORF"
    ## [1] "Analyzing TOP motif"
    ## [1] "Analyzing RNA structure"
    ## [1] "Export all data"

# Session information

``` r
sessionInfo()
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
    ##  [1] ORFik_1.8.1                       BSgenome.Hsapiens.UCSC.hg38_1.4.3
    ##  [3] BSgenome_1.56.0                   knitr_1.28                       
    ##  [5] stringr_1.4.0                     magrittr_1.5                     
    ##  [7] data.table_1.12.8                 dplyr_1.0.0                      
    ##  [9] khroma_1.3.0                      ggplot2_3.3.1                    
    ## [11] GenomicAlignments_1.24.0          Rsamtools_2.4.0                  
    ## [13] Biostrings_2.56.0                 XVector_0.28.0                   
    ## [15] SummarizedExperiment_1.18.1       DelayedArray_0.14.0              
    ## [17] matrixStats_0.56.0                GenomicFeatures_1.40.0           
    ## [19] AnnotationDbi_1.50.0              Biobase_2.48.0                   
    ## [21] rtracklayer_1.48.0                GenomicRanges_1.40.0             
    ## [23] GenomeInfoDb_1.24.0               IRanges_2.22.1                   
    ## [25] S4Vectors_0.26.0                  BiocGenerics_0.34.0              
    ## [27] rmarkdown_2.2                    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2             splines_4.0.0          bit64_0.9-7           
    ##  [4] assertthat_0.2.1       askpass_1.1            BiocFileCache_1.12.0  
    ##  [7] blob_1.2.1             GenomeInfoDbData_1.2.3 yaml_2.2.1            
    ## [10] progress_1.2.2         pillar_1.4.4           RSQLite_2.2.0         
    ## [13] lattice_0.20-41        glue_1.4.1             digest_0.6.25         
    ## [16] RColorBrewer_1.1-2     colorspace_1.4-1       plyr_1.8.6            
    ## [19] htmltools_0.4.0        Matrix_1.2-18          DESeq2_1.28.0         
    ## [22] XML_3.99-0.3           pkgconfig_2.0.3        biomaRt_2.44.0        
    ## [25] genefilter_1.70.0      zlibbioc_1.34.0        xtable_1.8-4          
    ## [28] purrr_0.3.4            scales_1.1.1           BiocParallel_1.22.0   
    ## [31] annotate_1.66.0        tibble_3.0.1           openssl_1.4.1         
    ## [34] generics_0.0.2         ellipsis_0.3.1         withr_2.4.1           
    ## [37] survival_3.1-12        crayon_1.3.4           memoise_1.1.0         
    ## [40] evaluate_0.14          GGally_2.0.0           tools_4.0.0           
    ## [43] prettyunits_1.1.1      hms_0.5.3              lifecycle_0.2.0       
    ## [46] locfit_1.5-9.4         munsell_0.5.0          compiler_4.0.0        
    ## [49] rlang_0.4.10           grid_4.0.0             RCurl_1.98-1.2        
    ## [52] rappdirs_0.3.1         bitops_1.0-6           gtable_0.3.0          
    ## [55] reshape_0.8.8          DBI_1.1.0              curl_4.3              
    ## [58] R6_2.4.1               gridExtra_2.3          bit_1.1-15.2          
    ## [61] stringi_1.4.6          Rcpp_1.0.4.6           geneplotter_1.66.0    
    ## [64] vctrs_0.3.1            dbplyr_1.4.4           tidyselect_1.1.0      
    ## [67] xfun_0.14
