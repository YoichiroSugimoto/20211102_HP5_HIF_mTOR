s4-2-1 Transcript assignment to TSS (RCC-4) version 5
================
Yoichiro Sugimoto
02 March, 2022

  - [Overview](#overview)
  - [Define variables](#define-variables)
  - [Select the data to define transcript and TSS
    assignment](#select-the-data-to-define-transcript-and-tss-assignment)
  - [Assign transcript to TSS](#assign-transcript-to-tss)
  - [Session information](#session-information)

# Overview

A transcript will be assigned to a TSS based on the 5’ end-Seq of total
mRNA data. This will be performed for RCC4 and 786-O cells
independently, and both the cell lines with and without VHL will be
considered.

5’ end-Seq reads starting from each TSS will be assembled using
`StringTie` software with the parameters, `--conservative -j 5`. To
guide the assembly, this analysis will consider both TSS defined before
(the median position will be used as the input) and transcript
annotation. The analysis will generate possible reconstructed 5’ ends of
source transcripts (assemblies) for each TSS with the 5′ end-Seq read
coverage score (calculated by multiplying the transcript per million
(TPM) estimate and length of the assembly reported by the software). For
each TSS, any assembly with a lower coverage score (equal or less than
80% of the highest score) will be filtered out, since these are likely
to have less data support or from poorly expressed isoforms.

The selected assemblies are plausible 5’ ends of expressed transcripts
from a TSS defined by data. To complement the remaining part of
transcript information, for each assembly, high homology annotated
transcripts those fall into the following categories will be searched:
(i). for an assembly that comprises multiple exons, annotated
transcripts containing all the exon junctions of the assembly except the
donor of the first exon junction (categorised as i-a when the donor of
the first exon junction is also contained, and otherwise categorised as
i-b); (ii) for an assembly comprised of one exon, annotated transcripts
containing 3’ end of the assembly.

Following this process, it remains possible that more than one
assembly-annotated transcript pairs are associated with a TSS. In this
case, we will prioritize the assembly-transcript pairs using the
following criteria:

1.  transcripts categorized as protein coding
2.  start codon of transcripts downstream to TSS (q50)
3.  prioritize assembly-transcript pairs by the overlap category
    (described above) in the following order: (i-a), (ii), and (i-b)
4.  prioritize transcripts annotated by RefSeq entries over GENCODE
5.  transcripts flagged as “basic” transcripts by GENCODE
6.  5’ UTR length of transcripts not equal to 0
7.  3’ UTR length of transcripts not equal to 0
8.  transcripts with long CDS length
9.  transcripts with long 5’ UTR
10. transcripts with long 3’ UTR, and
11. assemblies with high coverage score.

In this way, an annotated transcript will be assigned to a TSS.

Finally, TSS of the assigned transcripts will be corrected to the
experimentally defined TSS. If the TSS is downstream to the start codon
of the assigned transcripts, 5’ most AUG sequence with the same reading
frame as the original CDS will be searched and the region between the
AUG sequence and the original stop codon will be defined as CDS of the
mRNA. The corrected transcripts will be used as the assigned transcripts
to the TSS.

``` r
## Bioconductor packages
library("rtracklayer")
library("GenomicFeatures")
library("GenomicAlignments")
## Specify the number of CPUs to be used
processors <- 16

temp <- sapply(list.files("../functions", full.names = TRUE), source)
source("./functions/assignTxToTss.R")

sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

annot.dir <- normalizePath(file.path("../../annotation/"))
hg38.annot.dir <- file.path(
    annot.dir,
    "hg38_annotation"
)
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")

software.path.dt <- file.path("../../data/environment/software_path.csv") %>% fread

ref.seq.dir <- file.path(annot.dir, "hg38_annotation/ref_sequences/")
genome.fa.file <- list.files(ref.seq.dir, pattern = "_genome.fa$", full.names = TRUE)

annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)

kallisto.index.dir <- file.path(annot.dir, "hg38_annotation/kallisto_indices/")

results.dir <- file.path("../../results")
s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
star.aligned.bam.dir <- file.path(s2.alignment.dir, "s2-1-b-star-aligned_bam")
s2.2.processed.bam.dir <-  file.path(s2.alignment.dir, "s2-2-processed-data")
s2.2.1.tss.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-1-tss-bam")

## Previous results
s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")

s4.2.tx.assignment.dir <- file.path(s4.tss.dir, "s4-2-transcript-assignment")
s4.2.1.tss.tx.map.RCC4.dir <- file.path(s4.2.tx.assignment.dir, "s4-2-1-tss-transcript-mapping-RCC4")

create.dirs(c(
    s4.2.tx.assignment.dir,
    s4.2.1.tss.tx.map.RCC4.dir
))

sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]
```

# Define variables

``` r
filtered.tss.with.quantile.file <- file.path(
  s4.1.6.filtered.tss.dir,
  "filtered-tss-with-quantile.csv"
)

filtered.tss.with.quantile.dt <- fread(filtered.tss.with.quantile.file)


all.tx.gtf <- file.path(
    annot.ps.dir,
    "all-transcript.gtf"
)

all.primary.tx.info.file <- file.path(
    annot.ps.dir,
    "all_GENCODE_RefSeq_transcript_info.csv"
)

all.primary.tx.dt <- fread(all.primary.tx.info.file)

tx.map.out.dir <- s4.2.1.tss.tx.map.RCC4.dir

function.dir <- file.path("./functions")
```

# Select the data to define transcript and TSS assignment

``` r
## Below is for RCC4 data
sl.sample.names <- grep("total", sample.names, value = TRUE) %>%
  grep("_HIF1B_N_", ., value = TRUE) %>%
  grep("RCC4_noVHL|RCC4_VHL", ., value = TRUE)

print("The following data will be used to assign transcripts to TSS")
```

    ## [1] "The following data will be used to assign transcripts to TSS"

``` r
print(sl.sample.names)
```

    ## [1] "total_RCC4_VHL_HIF1B_N_1"   "total_RCC4_VHL_HIF1B_N_3"  
    ## [3] "total_RCC4_VHL_HIF1B_N_4"   "total_RCC4_noVHL_HIF1B_N_1"
    ## [5] "total_RCC4_noVHL_HIF1B_N_3" "total_RCC4_noVHL_HIF1B_N_4"

# Assign transcript to TSS

``` r
temp <- assignTxToTss(
    tx.map.out.dir = tx.map.out.dir,
    filtered.tss.with.quantile.dt = filtered.tss.with.quantile.dt,
    all.tx.gtf = all.tx.gtf,
    all.primary.tx.dt = all.primary.tx.dt,
    sl.sample.names = sl.sample.names,
    star.aligned.bam.dir = star.aligned.bam.dir,
    software.path.dt = software.path.dt,
    cell.name = "RCC4",
    function.dir = function.dir
)
```

    ## [1] "The following sample groups are examined:"
    ## [1] "RCC4_VHL_HIF1B_N"   "RCC4_noVHL_HIF1B_N"
    ## [1] "Bam files are split by TSS"
    ## [1] "Run StringTie"
    ## [1] "Assign Tx to TSS"
    ## [1] "Export GTF for assigned Tx"

    ## Loading required package: BSgenome

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    ## [1] "Exporting data"

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
