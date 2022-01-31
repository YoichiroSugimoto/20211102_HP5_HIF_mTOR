---
title: "s6-1 Gene level mRNA differential expression analysis"
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

Gene level mRNA differential expression analysis will be performed here.


```r
## Biocoductor packages
library("DESeq2")
## Data analysis packages
library("magrittr")
library("matrixStats")

## Parallelization
## Specify the number of CPUs to be used
processors <- 8
library("BiocParallel")
register(MulticoreParam(processors))

temp <- sapply(list.files("../functions", full.names = TRUE), source)
source(file.path("./functions/load_total_sample_data.R"), chdir = TRUE)
```

```
## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
## [1] "The following R objects were exported: total.sample.dt, total.coldata.df, total.comparison.dt"
```

```r
source(file.path("./functions/test_differential_expression.R"), chdir = TRUE)
```

```
## [1] "The following functions are exported: calcMeanCount(), analyzeDE()"
```

```r
set.seed(0)
```



```r
annot.dir <- normalizePath(file.path("../../annotation/"))
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)

results.dir <- file.path("../../results")
s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")
s4.1.7.count.per.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-7-count-per-tss")
s4.2.tx.assignment.dir <- file.path(s4.tss.dir, "s4-2-transcript-assignment")
s4.2.1.tss.tx.map.RCC4.dir <- file.path(s4.2.tx.assignment.dir, "s4-2-1-tss-transcript-mapping-RCC4")

s6.dir <- file.path(results.dir, "s6-differential-regulation-analysis")
s6.1.dir <- file.path(s6.dir, "s6-1-differentially-expressed-genes")

create.dirs(c(
  s6.dir,
  s6.1.dir
))
```

# Data preparation


```r
filtered.tss.with.quantile.dt <- file.path(
  s4.1.6.filtered.tss.dir,
  "filtered-tss-with-quantile.csv"
) %>% fread

total.count.per.tss.dt <- file.path(s4.1.7.count.per.tss.dir, "count-per-confident-tss.csv") %>%
    fread

total.count.per.tss.dt <- total.count.per.tss.dt[
  , c("tss_name", grep("^total_", colnames(total.count.per.tss.dt), value = TRUE)),
    with = FALSE
]

total.count.dt <- countByGeneFromTss(total.count.per.tss.dt)

colnames(total.count.dt) <- gsub(
    "total_", "", colnames(total.count.dt)
)

total.count.dt <- total.count.dt[, c("gene_id", rownames(total.coldata.df)), with = FALSE]

total.count.dt <- merge(
    primary.tx.dt[!duplicated(gene_id), .(gene_id, gene_name, biotype)],
    total.count.dt,
    by = "gene_id"
)
```

# Analysis of differentially expressed genes

## Data preparation for DESeq2 analysis


```r
countDt2Df <- function(count.dt){
    ## temp.count.dt <- copy(count.dt[biotype == "protein_coding"])
    temp.count.dt <- copy(count.dt[biotype != ""])
    count.df <- as.data.frame(temp.count.dt)
    rownames(count.df) <- count.df[, "gene_id"]
    annot.df <- count.df[, c("gene_id", "gene_name", "biotype")]
    count.df <- dropColumnDf(count.df, drop.vec = c("gene_id", "gene_name", "biotype"))
    return(list(count.df = count.df, annot.df = annot.df))
}

deseq2.in.list <- countDt2Df(total.count.dt)

if(!all(colnames(deseq2.in.list$count.df) == rownames(total.coldata.df))){
    stop("colnames of count.df does not match with rownames of coldata.df")
} else {"OK"}
```

```
## [1] "OK"
```

## List of comparison


```r
## This is now performed at the R file in functions/ directory
print("List of comparisons")
```

```
## [1] "List of comparisons"
```

```r
total.comparison.dt
```

```
##                         comparison     exp_formula
## 1:   RCC4_xx_HIF1B_N__noVHL_vs_VHL            ~VHL
## 2:   786O_xx_HIF1B_N__noVHL_vs_VHL            ~VHL
## 3: 786O_xx_noHIF1B_N__noVHL_vs_VHL            ~VHL
## 4:     RCC4_noVHL_HIF1B_xx__H_vs_N ~clone + oxygen
## 5:       RCC4_VHL_HIF1B_xx__H_vs_N ~clone + oxygen
```


## Analysis of differential gene expression by VHL or oxygen


```r
de.res.dts <- mapply(
    analyzeDE,
    total.count.dt = list(total.count.dt),
    annot.dt = list(data.table(deseq2.in.list$annot.df)),
    ref.column.name = list("gene_id"),
    input.sample.data.df = list(total.coldata.df),
    comparison.name = total.comparison.dt[, comparison],
    exp.design = total.comparison.dt[, exp_formula],
    out.dir = list(s6.1.dir),
    SIMPLIFY = FALSE
)
```

```
## [1] "Analyzing RCC4_xx_HIF1B_N__noVHL_vs_VHL"
## [1] "With the formula of ~VHL"
## DataFrame with 6 rows and 6 columns
##                          cell      VHL    HIF1B   oxygen    clone sizeFactor
##                      <factor> <factor> <factor> <factor> <factor>  <numeric>
## RCC4_VHL_HIF1B_N_1       RCC4    VHL      HIF1B        N       11   1.106084
## RCC4_VHL_HIF1B_N_3       RCC4    VHL      HIF1B        N       13   1.079053
## RCC4_VHL_HIF1B_N_4       RCC4    VHL      HIF1B        N       14   1.076944
## RCC4_noVHL_HIF1B_N_1     RCC4    noVHL    HIF1B        N       1    1.203965
## RCC4_noVHL_HIF1B_N_3     RCC4    noVHL    HIF1B        N       3    0.874438
## RCC4_noVHL_HIF1B_N_4     RCC4    noVHL    HIF1B        N       4    0.763353
```

![](s6-1-gene-level-mRNA-abundance-change-analysis_files/figure-html/Analsyis of VHL dependent mRNA abundance changes-1.png)<!-- -->

```
## [1] "Summary of analysis results"
##                  reg_direction
## stat_significance    Up  Down   Sum
##   padj < 0.1       1014   987  2001
##   not significant  5892  6075 11967
##   Sum              6906  7062 13968
## [1] "Analyzing 786O_xx_HIF1B_N__noVHL_vs_VHL"
## [1] "With the formula of ~VHL"
## DataFrame with 8 rows and 6 columns
##                          cell      VHL    HIF1B   oxygen    clone sizeFactor
##                      <factor> <factor> <factor> <factor> <factor>  <numeric>
## 786O_VHL_HIF1B_N_1       786O    VHL      HIF1B        N       11   1.065814
## 786O_VHL_HIF1B_N_2       786O    VHL      HIF1B        N       12   1.107852
## 786O_VHL_HIF1B_N_3       786O    VHL      HIF1B        N       13   0.999558
## 786O_VHL_HIF1B_N_4       786O    VHL      HIF1B        N       14   1.261652
## 786O_noVHL_HIF1B_N_1     786O    noVHL    HIF1B        N       1    1.040897
## 786O_noVHL_HIF1B_N_2     786O    noVHL    HIF1B        N       2    0.988377
## 786O_noVHL_HIF1B_N_3     786O    noVHL    HIF1B        N       3    0.841843
## 786O_noVHL_HIF1B_N_4     786O    noVHL    HIF1B        N       4    0.836513
```

![](s6-1-gene-level-mRNA-abundance-change-analysis_files/figure-html/Analsyis of VHL dependent mRNA abundance changes-2.png)<!-- -->

```
## [1] "Summary of analysis results"
##                  reg_direction
## stat_significance    Up  Down   Sum
##   padj < 0.1        137   207   344
##   not significant  6240  6065 12305
##   Sum              6377  6272 12649
## [1] "Analyzing 786O_xx_noHIF1B_N__noVHL_vs_VHL"
## [1] "With the formula of ~VHL"
## DataFrame with 6 rows and 6 columns
##                            cell      VHL    HIF1B   oxygen    clone sizeFactor
##                        <factor> <factor> <factor> <factor> <factor>  <numeric>
## 786O_VHL_noHIF1B_N_1       786O    VHL    noHIF1B        N      111   0.932247
## 786O_VHL_noHIF1B_N_2       786O    VHL    noHIF1B        N      112   0.854477
## 786O_VHL_noHIF1B_N_3       786O    VHL    noHIF1B        N      113   0.708944
## 786O_noVHL_noHIF1B_N_1     786O    noVHL  noHIF1B        N      101   1.352837
## 786O_noVHL_noHIF1B_N_2     786O    noVHL  noHIF1B        N      102   1.114647
## 786O_noVHL_noHIF1B_N_3     786O    noVHL  noHIF1B        N      103   1.245521
```

```
## [1] "Summary of analysis results"
##                  reg_direction
## stat_significance    Up  Down   Sum
##   padj < 0.1          4     1     5
##   not significant  7178  6919 14097
##   Sum              7182  6920 14102
```

```
## factor levels were dropped which had no samples
```

![](s6-1-gene-level-mRNA-abundance-change-analysis_files/figure-html/Analsyis of VHL dependent mRNA abundance changes-3.png)<!-- -->

```
## [1] "Analyzing RCC4_noVHL_HIF1B_xx__H_vs_N"
## [1] "With the formula of ~clone + oxygen"
## DataFrame with 6 rows and 6 columns
##                          cell      VHL    HIF1B   oxygen    clone sizeFactor
##                      <factor> <factor> <factor> <factor> <factor>  <numeric>
## RCC4_noVHL_HIF1B_N_1     RCC4    noVHL    HIF1B        N        1   1.201969
## RCC4_noVHL_HIF1B_N_3     RCC4    noVHL    HIF1B        N        3   0.872924
## RCC4_noVHL_HIF1B_N_4     RCC4    noVHL    HIF1B        N        4   0.761090
## RCC4_noVHL_HIF1B_H_1     RCC4    noVHL    HIF1B        H        1   1.386723
## RCC4_noVHL_HIF1B_H_3     RCC4    noVHL    HIF1B        H        3   0.955275
## RCC4_noVHL_HIF1B_H_4     RCC4    noVHL    HIF1B        H        4   0.966994
```

```
## [1] "Summary of analysis results"
##                  reg_direction
## stat_significance   Up Down  Sum
##   padj < 0.1       158  116  274
##   not significant 2239 1792 4031
##   Sum             2397 1908 4305
```

```
## factor levels were dropped which had no samples
```

![](s6-1-gene-level-mRNA-abundance-change-analysis_files/figure-html/Analsyis of VHL dependent mRNA abundance changes-4.png)<!-- -->

```
## [1] "Analyzing RCC4_VHL_HIF1B_xx__H_vs_N"
## [1] "With the formula of ~clone + oxygen"
## DataFrame with 6 rows and 6 columns
##                        cell      VHL    HIF1B   oxygen    clone sizeFactor
##                    <factor> <factor> <factor> <factor> <factor>  <numeric>
## RCC4_VHL_HIF1B_N_1     RCC4      VHL    HIF1B        N       11   0.961696
## RCC4_VHL_HIF1B_N_3     RCC4      VHL    HIF1B        N       13   0.939595
## RCC4_VHL_HIF1B_N_4     RCC4      VHL    HIF1B        N       14   0.940142
## RCC4_VHL_HIF1B_H_1     RCC4      VHL    HIF1B        H       11   1.224317
## RCC4_VHL_HIF1B_H_3     RCC4      VHL    HIF1B        H       13   1.142949
## RCC4_VHL_HIF1B_H_4     RCC4      VHL    HIF1B        H       14   0.867809
```

![](s6-1-gene-level-mRNA-abundance-change-analysis_files/figure-html/Analsyis of VHL dependent mRNA abundance changes-5.png)<!-- -->

```
## [1] "Summary of analysis results"
##                  reg_direction
## stat_significance    Up  Down   Sum
##   padj < 0.1       1157   775  1932
##   not significant  4142  4582  8724
##   Sum              5299  5357 10656
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
##  [1] knitr_1.28                  stringr_1.4.0              
##  [3] data.table_1.12.8           dplyr_1.0.0                
##  [5] khroma_1.3.0                ggplot2_3.3.1              
##  [7] BiocParallel_1.22.0         magrittr_1.5               
##  [9] DESeq2_1.28.0               SummarizedExperiment_1.18.1
## [11] DelayedArray_0.14.0         matrixStats_0.56.0         
## [13] Biobase_2.48.0              GenomicRanges_1.40.0       
## [15] GenomeInfoDb_1.24.0         IRanges_2.22.1             
## [17] S4Vectors_0.26.0            BiocGenerics_0.34.0        
## [19] rmarkdown_2.2              
## 
## loaded via a namespace (and not attached):
##  [1] bdsmatrix_1.3-4        Rcpp_1.0.4.6           locfit_1.5-9.4        
##  [4] mvtnorm_1.1-1          apeglm_1.10.0          lattice_0.20-41       
##  [7] digest_0.6.25          plyr_1.8.6             R6_2.4.1              
## [10] emdbook_1.3.12         coda_0.19-3            RSQLite_2.2.0         
## [13] evaluate_0.14          pillar_1.4.4           zlibbioc_1.34.0       
## [16] rlang_0.4.10           annotate_1.66.0        blob_1.2.1            
## [19] Matrix_1.2-18          bbmle_1.0.23.1         splines_4.0.0         
## [22] geneplotter_1.66.0     RCurl_1.98-1.2         bit_1.1-15.2          
## [25] munsell_0.5.0          numDeriv_2016.8-1.1    compiler_4.0.0        
## [28] xfun_0.14              pkgconfig_2.0.3        htmltools_0.4.0       
## [31] tidyselect_1.1.0       tibble_3.0.1           GenomeInfoDbData_1.2.3
## [34] XML_3.99-0.3           crayon_1.3.4           withr_2.4.1           
## [37] MASS_7.3-51.6          bitops_1.0-6           grid_4.0.0            
## [40] xtable_1.8-4           gtable_0.3.0           lifecycle_0.2.0       
## [43] DBI_1.1.0              scales_1.1.1           stringi_1.4.6         
## [46] XVector_0.28.0         genefilter_1.70.0      ellipsis_0.3.1        
## [49] vctrs_0.3.1            generics_0.0.2         RColorBrewer_1.1-2    
## [52] tools_4.0.0            bit64_0.9-7            glue_1.4.1            
## [55] purrr_0.3.4            survival_3.1-12        yaml_2.2.1            
## [58] AnnotationDbi_1.50.0   colorspace_1.4-1       memoise_1.1.0
```
