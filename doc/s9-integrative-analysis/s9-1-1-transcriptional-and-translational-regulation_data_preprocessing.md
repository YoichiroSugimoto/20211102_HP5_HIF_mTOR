---
title: "s9-1-1 Intersection of transcriptional and translational regulation (data preprocessing)"
author: "Yoichiro Sugimoto"
date: "06 January, 2022"
vignette: >
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



```r
## Specify the number of CPUs to be used
processors <- 8

temp <- sapply(list.files("../functions", full.names = TRUE), source, chdir = TRUE)
temp <- sapply(list.files("../s8-analysis-of-translation/functions", full.names = TRUE), source, chdir = TRUE)
```

```
## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
## [1] "The following objects are exported: poly.coldata.df, poly.sample.dt, translation.comparison.dt"
## [1] "In translation.comparison.dt, xx specifies the factor compared where the comparison is specified after __, while yy is a wildcard. From left, each factor specifies cell, VHL, EIF4E2, clone, and treatment"
## [1] "The following functions were exported: analyzeDtg(), subsetColdata()"
```

```r
source("../s6-differential-expression-and-tss-usage/functions/load_total_analysis_results.R", chdir = TRUE)
```

```
## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
## [1] "The following R objects were exported: total.sample.dt, total.coldata.df, total.comparison.dt"
## [1] "Comparison information was loaded"
## [1] "/camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/results"
## [1] "The following objects were loaded: tss.de.res.dt, tss.ratio.res.dt, diff.tss.res.dt"
```

```r
set.seed(0)
```



```r
annot.dir <- normalizePath(file.path("../../annotation/"))
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
all.primary.tx.dt <- file.path(
    annot.ps.dir, "all_GENCODE_RefSeq_transcript_info.csv"
) %>% fread

results.dir <- file.path("../../results")

s6.dir <- file.path(results.dir, "s6-differential-regulation-analysis")
s6.1.dir <- file.path(s6.dir, "s6-1-differentially-expressed-genes")

s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs")
s8.1.1.dir <- file.path(s8.1.dir, "gene-level-dte")
s8.3.dir <- file.path(s8.dir, "s8-3-validation-of-method")

s9.dir <- file.path(results.dir, "s9-integrative-analysis")

create.dirs(
    c(
        s9.dir
    )
)

sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")
sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]
```

# Summarization of differential gene expression regulation analysis results

## Differential expression data (gene level) import



```r
## Differential gene expression regulation
readDeResults <- function(comparison.name, result.file.dir, analysis.level = "gene", log2fc_threshold = log2(1.5)){
    
    de.res.dt <- fread(file.path(
        result.file.dir,
        paste0(comparison.name, "_DE.csv")
    ))
    de.res.dt[, `:=`(
        mRNA_regulation = case_when(
            padj < 0.1 & log2fc > log2fc_threshold ~ "Up",
            padj < 0.1 & log2fc < -log2fc_threshold ~ "Down",
            TRUE ~ "Not_significant" 
        ),
        analysis_level = analysis.level,
        comparison_name = comparison.name
    )]

    data.cols <- c(
        "log2fc", "shrlog2fc", "padj",
        "meanNormCount_treated", "meanNormCount_base",
        "mRNA_regulation"
    )

    de.res.dt <- de.res.dt[, c("gene_id", data.cols), with = FALSE]

    setnames(
        de.res.dt,
        old = data.cols,
        new = paste0(data.cols, "_", comparison.name)
    )
    
    return(de.res.dt)
}

de.gene.res.dts <- lapply(
    total.comparison.dt[, comparison],
    readDeResults,
    result.file.dir = s6.1.dir,
    analysis.level = "gene"
) 

de.gene.res.dt <- Reduce(
    function(...) merge(..., all = TRUE, by = "gene_id"), de.gene.res.dts
)
```

## Differential translation data import


```r
readTranslationAnalysisRes <- function(dte.comparison.name, result.file.dir, analysis.level = "gene"){
    dte.comparison.filename <- dte.comparison.name %>%
        gsub("\\(", "", .) %>%
        gsub("\\)", "", .) %>%
        gsub("\\|", "-", .)
    
    dte.res.dt <- fread(file.path(
        result.file.dir,
        paste0(dte.comparison.filename, ".csv")
    ))
    dte.res.dt[, `:=`(
        analysis_level = analysis.level,
        comparison_name = dte.comparison.name
    )]
    
    return(dte.res.dt)    
}

all.trsl.res.dt <- lapply(
    translation.comparison.dt[1:6, comparison],
    readTranslationAnalysisRes,
    result.file.dir = s8.1.1.dir,
    analysis.level = "gene"
) %>% rbindlist(use.names = TRUE, fill = TRUE)

d.all.trsl.res.dt <- dcast(
    all.trsl.res.dt,
    gene_id ~ comparison_name,
    value.var = c(
        "MRL_log2fc",
        "padj_translation",
        "MRL_treated", "MRL_base", "translational_regulation"
    )
)

## Merge all the analysis results
all.de.dte.res.dt <- merge(
    all.primary.tx.dt[!duplicated(gene_id), .(gene_id, gene_name, biotype, chromosome_name)],
    de.gene.res.dt,
    by = "gene_id",
    all.y = TRUE
) %>%
    merge(
        d.all.trsl.res.dt,
        by = "gene_id",
        all = TRUE
    )

all.de.dte.res.dt <- all.de.dte.res.dt[biotype == "protein_coding"]
```


# Classifications of genes by the VHL and mTOR dependent regulation



```r
all.de.dte.res.dt[, `:=`(
    VHL_target_RCC4 = case_when(
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Up" ~ "VHL_loss_induced",
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Down" ~ "VHL_loss_repressed",
        TRUE ~ "non_VHL_target"
    ) %>% factor(levels = c(
                     "VHL_loss_repressed",
                     "non_VHL_target",
                     "VHL_loss_induced"
                 )),
    RCC4_mRNA_mTOR_trsl_group = case_when(
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Up" &
        translational_regulation_RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA == "Up" ~
            "mRNA_up_and_translation_up",
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Up" &
        translational_regulation_RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA == "Down" ~
            "mRNA_up_and_translation_down",
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Down" &
        translational_regulation_RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA == "Up" ~
            "mRNA_down_and_translation_up",
        mRNA_regulation_RCC4_xx_HIF1B_N__noVHL_vs_VHL == "Down" &
        translational_regulation_RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA == "Down" ~
            "mRNA_down_and_translation_down",
        TRUE ~ "Others"
    )
)]

all.de.dte.res.dt[, table(VHL_target_RCC4)]
```

```
## VHL_target_RCC4
## VHL_loss_repressed     non_VHL_target   VHL_loss_induced 
##                508              10577                433
```

```r
all.de.dte.res.dt[, table(RCC4_mRNA_mTOR_trsl_group)]
```

```
## RCC4_mRNA_mTOR_trsl_group
## mRNA_down_and_translation_down   mRNA_down_and_translation_up 
##                            103                             77 
##   mRNA_up_and_translation_down     mRNA_up_and_translation_up 
##                            124                            131 
##                         Others 
##                          11083
```


# Annotation of the functions of genes



```r
key.go.gene.dt <- file.path(
    s8.3.dir,
    "key_GO_genes.csv"
) %>%
    fread

all.de.dte.res.dt[, `:=`(
    mRNA_trsl_intersection_by_functions = case_when(
        gene_id %in% key.go.gene.dt[GO_name == "Glycolysis", gene_id] ~ "Glycolysis",
        gene_id %in% key.go.gene.dt[
                         GO_name %in%
                         c("Vascular Process", "Angiogenesis"), gene_id] ~
            "Angiogenesis or vascular process",
        TRUE ~ "Others"
    ) %>% factor(
              levels = c("Glycolysis", "Angiogenesis or vascular process", "Others")
          )
)]

all.de.dte.res.dt[, table(mRNA_trsl_intersection_by_functions)]
```

```
## mRNA_trsl_intersection_by_functions
##                       Glycolysis Angiogenesis or vascular process 
##                               42                              371 
##                           Others 
##                            11105
```

```r
fwrite(
    all.de.dte.res.dt,
    file.path(
        s9.dir,
        "all-differential-expression-and-translation-data.csv"
    )
)
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
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] knitr_1.28        stringr_1.4.0     magrittr_1.5      data.table_1.12.8
## [5] dplyr_1.0.0       khroma_1.3.0      ggplot2_3.3.1     rmarkdown_2.2    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6     tidyselect_1.1.0 munsell_0.5.0    colorspace_1.4-1
##  [5] R6_2.4.1         rlang_0.4.10     tools_4.0.0      grid_4.0.0      
##  [9] gtable_0.3.0     xfun_0.14        withr_2.4.1      htmltools_0.4.0 
## [13] ellipsis_0.3.1   yaml_2.2.1       digest_0.6.25    tibble_3.0.1    
## [17] lifecycle_0.2.0  crayon_1.3.4     purrr_0.3.4      vctrs_0.3.1     
## [21] glue_1.4.1       evaluate_0.14    stringi_1.4.6    compiler_4.0.0  
## [25] pillar_1.4.4     generics_0.0.2   scales_1.1.1     pkgconfig_2.0.3
```
