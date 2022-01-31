---
title: "s8-3-2-1 Analysis of mTOR-dependent translational regulation (1/2)"
author: "Yoichiro Sugimoto"
date: "30 January, 2022"
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


The translational effect of mTOR inhibition in RCC-4 VHL will be examined here.



```r
## Specify the number of CPUs to be used
processors <- 8
## library("BiocParallel")
## register(MulticoreParam(processors))

temp <- sapply(list.files("../functions", full.names = TRUE), source)
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
sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

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
s4.3.tx.info.dir <- file.path(s4.tss.dir, "s4-3-transcript-info")
s4.3.1.tx.info.rcc4.dir <- file.path(s4.3.tx.info.dir, "s4-3-1-transcript-info-for-RCC4")

s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs")
s8.1.1.dir <- file.path(s8.1.dir, "gene-level-dte")
s8.1.2.dir <- file.path(s8.1.dir, "tx-level-dte")
s8.3.dir <- file.path(s8.dir, "s8-3-validation-of-method")

sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]
```

# Evaluation of translation changes by mTOR inhibition (gene level)

## VHL status for the analysis



```r
vhl.status <- "VHL"
```


## Data import


```r
all.filtered.gene.dt <- file.path(
    s8.3.dir,
    "filtered_gene_for_polysome_analysis.csv"
) %>% fread


all.filtered.gene.dt[RCC4_VHL_NA == TRUE & RCC4_VHL_Torin1 == TRUE][
    gene_id %in%
    primary.tx.dt[!duplicated(gene_id)][biotype == "protein_coding", gene_id]    
] %>%
    nrow
```

```
## [1] 8952
```




```r
torin.gene.trsl.dt <- file.path(
    s8.1.1.dir,
    paste0(
        "RCC4_",
        vhl.status,
        "_EIF4E2_yy_xx__Torin1_vs_NA.csv"
    )
) %>% fread

torin.gene.trsl.dt[
  , trsl_reg_class := case_when(
        translational_regulation %in% "Up" ~ "Preserved",
        translational_regulation %in% "Down" ~ "Repressed",
        TRUE ~ "Not significant"
    ) %>%
        factor(levels = c("Preserved", "Not significant", "Repressed"))
]
```

## Analysis of the absolute changes in translation upon mTOR inhibition



```r
torin.gene.trsl.dt <- torin.gene.trsl.dt[
    gene_id %in% all.filtered.gene.dt[RCC4_VHL_NA == TRUE & RCC4_VHL_Torin1 == TRUE, gene_id] &
    gene_id %in% primary.tx.dt[!duplicated(gene_id)][biotype == "protein_coding", gene_id]
]

ggplot(
    data = torin.gene.trsl.dt,
    aes(
        x = MRL_base,
        y = MRL_treated
    )
) +
    geom_abline(slope = 1, intercept = 0, color = "gray60") +
    geom_point() +
    xlab("MRL without Torin 1") +
    ylab("MRL with Torin 1") +
    xlim(c(1, 8)) + ylim(c(1, 8)) +
    theme(
        aspect.ratio = 1,
        legend.position = "bottom"
    )
```

![](s8-3-2-1-analysis-of-mTOR-dependent-translational-regulation-1_files/figure-html/plot differential translation data-1.png)<!-- -->


# Analysis of changes in translational efficiency as a function of the functional classes of mRNAs



```r
## Plot translation changes by mRNA classes
plotTrslChangeByClass <- function(trsl.dt, s8.3.dir, plot.ylab, sig.th = 0.05){
    ## trsl.dt must have 2 columns: gene_id and MRL_log2fc

    trsl.dt <- trsl.dt[!is.na(MRL_log2fc)]

    kegg.all.dt <- file.path(s8.3.dir, "key_KEGG_genes.csv") %>%
        fread

    all.sl.terms <- kegg.all.dt[, unique(term_name)]
    all.sl.term.groups <- kegg.all.dt[, unique(term_group)]

    mrna.class.trsl.dt <- merge(
        kegg.all.dt,
        trsl.dt,
        by = "gene_id"
    )

    all.trsl.dt <- copy(trsl.dt) %>%
        {.[, `:=`(
             term_name = "All",
             term_group = "All"
         )]}

    mrna.class.trsl.dt <- rbind(
        mrna.class.trsl.dt,
        all.trsl.dt,
        use.names = TRUE
    )

    calcWilP <- function(sl.term, mrna.class.trsl.dt){
        dt <- data.table(
            term_name = sl.term,
            p_value = wilcox.test(
                mrna.class.trsl.dt[term_name == sl.term, MRL_log2fc],
                mrna.class.trsl.dt[
                    term_name == "All" &
                    !(gene_id %in% mrna.class.trsl.dt[term_name == sl.term, gene_id]),
                    MRL_log2fc
                ],
                alternative = "two.sided"
            )$p.value,
            rg =
                mrna.class.trsl.dt[order(
                    term_name == "All"
                )][term_name %in% c(sl.term, "All")][
                    !duplicated(gene_id)
                ] %$%
                rcompanion::wilcoxonRG(
                                x = MRL_log2fc,
                                g = term_name == "All"
                            ),
            N = nrow(mrna.class.trsl.dt[term_name == sl.term]),
            N_all = nrow(mrna.class.trsl.dt[term_name == "All"])
        )
        
        return(dt)
    }

    all.term.sig.dt <- lapply(
        all.sl.terms,
        calcWilP,
        mrna.class.trsl.dt = mrna.class.trsl.dt
    ) %>%
        rbindlist

    all.term.sig.dt[, `:=`(
        padj = p.adjust(p_value, method = "holm")
    )]

    all.term.sig.dt[, `:=`(
        direction = case_when(
            padj < sig.th & rg > 0 ~ "Up",
            padj < sig.th & rg < 0 ~ "Down",
            TRUE ~ "N.S."
        ),
        sig_mark = case_when(
            padj < sig.th * 0.1 ~ "**",
            padj < sig.th ~ "*",
            TRUE ~ NA_character_
        )
    )]
    
    print(all.term.sig.dt)

    all.term.sig.dt <- rbind(
        all.term.sig.dt,
        data.table(
            term_name = "All", direction = "All",
            N = nrow(mrna.class.trsl.dt[term_name == "All"])
        ),
        use.names = TRUE, fill = TRUE
    )

    all.term.sig.dt[, term_name_n := paste0(
                          term_name, " (", N, ")"
                      )]
    
    mrna.class.trsl.dt <- merge(
        mrna.class.trsl.dt,
        all.term.sig.dt,
        by = "term_name"
    )

    mrna.class.trsl.dt[, `:=`(
        term_name = factor(term_name, levels = c(all.sl.terms, "All")),
        term_group = factor(term_group, levels = c(all.sl.term.groups, "All")),
        term_name_n = factor(term_name_n, levels = all.term.sig.dt[, term_name_n])
    )]
    
    g1 <- ggplot(
        data = mrna.class.trsl.dt,
        aes(
            x = term_name_n,
            y = MRL_log2fc,
            color = direction,
            fill = direction
        )
    ) +
        geom_hline(
            yintercept = median(trsl.dt[, MRL_log2fc], na.rm = TRUE)
        ) +
        geom_boxplot(outlier.shape = NA) +
        stat_summary(
            geom = 'text', aes(label = sig_mark),
            fun = function(x){boxplot.stats(x)$stats[5]}, 
            vjust = -0.8, color = "black", size = 5
        ) +
        facet_grid(
            ~ term_group,
            scales = "free",
            space = "free"
        ) +
        scale_color_manual(
            values = c(
                "Up" = "#4477AA", "Down" = "#EE6677",
                "N.S." = "gray40", "All" = "black"
            )
        ) +
        scale_fill_manual(
            values = c(
                "Up" = "lightsteelblue2", "Down" = "mistyrose",
                "N.S." = "white", "All" = "white"
            )
        ) +
        xlab("Functional class") +
        ylab(plot.ylab) +
        theme(
            legend.position = "bottom"
        ) +
        scale_x_discrete(guide = guide_axis(angle = 90))

    print(g1)

    return(all.term.sig.dt)
}


all.term.sig.dt <- plotTrslChangeByClass(
    trsl.dt = torin.gene.trsl.dt,
    s8.3.dir = s8.3.dir,
    plot.ylab = "MRL log2 fold change with Torin 1",
    sig.th = 0.05
)
```

```
##                            term_name      p_value      rg   N N_all
##  1:            Transcription factors 2.121055e-33  0.3240 486  8952
##  2:          Transcription machinery 5.725644e-01  0.0236 194  8952
##  3:         Messenger RNA biogenesis 8.268102e-01  0.0071 329  8952
##  4:                      Spliceosome 4.293154e-01 -0.0265 307  8952
##  5:             Cytoplasmic ribosome 8.007361e-32 -0.7860  75  8952
##  6:           Mitochondrial ribosome 1.379340e-13 -0.4990  74  8952
##  7:              Translation factors 1.142642e-06 -0.3140  81  8952
##  8: Chaperones and folding catalysts 7.632713e-09 -0.2610 166  8952
##  9:             Membrane trafficking 6.310913e-02 -0.0378 898  8952
## 10:                 Ubiquitin system 7.086448e-08  0.1400 526  8952
## 11:                       Proteasome 7.389658e-13 -0.5360  60  8952
## 12:                       Glycolysis 1.486860e-08 -0.5250  39  8952
## 13:        Pentose phosphate pathway 4.652497e-08 -0.7060  20  8952
## 14:                        TCA cycle 2.192175e-04 -0.4360  24  8952
## 15:          Fatty acid biosynthesis 1.822160e-03 -0.3410  28  8952
## 16:           Fatty acid degradation 2.470674e-03 -0.3200  30  8952
## 17:                           Oxphos 2.451842e-07 -0.2970 102  8952
## 18:            Nucleotide metabolism 1.138032e-05 -0.2780  84  8952
## 19:            Amino acid metabolism 1.143919e-06 -0.2220 163  8952
##             padj direction sig_mark
##  1: 4.030005e-32        Up       **
##  2: 1.000000e+00      N.S.     <NA>
##  3: 1.000000e+00      N.S.     <NA>
##  4: 1.000000e+00      N.S.     <NA>
##  5: 1.441325e-30      Down       **
##  6: 2.344879e-12      Down       **
##  7: 1.142642e-05      Down       **
##  8: 1.144907e-07      Down       **
##  9: 2.524365e-01      N.S.     <NA>
## 10: 8.503738e-07        Up       **
## 11: 1.182345e-11      Down       **
## 12: 2.081604e-07      Down       **
## 13: 6.048246e-07      Down       **
## 14: 1.534523e-03      Down       **
## 15: 1.093296e-02      Down        *
## 16: 1.235337e-02      Down        *
## 17: 2.697026e-06      Down       **
## 18: 9.104257e-05      Down       **
## 19: 1.142642e-05      Down       **
```

```
## Warning: Removed 5 rows containing missing values (geom_text).
```

![](s8-3-2-1-analysis-of-mTOR-dependent-translational-regulation-1_files/figure-html/torin by class-1.png)<!-- -->

# Translational targets of mTOR defined by previous studies


Known mTOR hypersensitive classes were defined based on previous genome-wide studies.



```r
known.mtor.target.dt <- file.path(
    "../../data/others/20201127_previously_reported_mTOR_target_genes.csv"
) %>%
    fread

known.mtor.target.dt <- known.mtor.target.dt[
    (method == "Ribosome profiling"  & reported_gene_name != "C3ORF38") |
    (method == "Polysome profiling" & PP242_log2fc < 0 & PP242_FDR < 0.05)
]

kegg.all.dt <- file.path(s8.3.dir, "key_KEGG_genes.csv") %>%
    fread

kegg.all.dt[, term_count := .N, by = term_name]

known.mtor.target.kegg.dt <- merge(
    kegg.all.dt,
    known.mtor.target.dt[, .(gene_id, gene_name, method)],
    by = "gene_id"
)

known.mtor.target.kegg.dt[, target_count := .N, by = list(term_name, method)]

target.count.per.term.dt <- known.mtor.target.kegg.dt[
    !duplicated(term_name, method),
    .(term_group, term_name, method, target_count, term_count)
]

target.count.per.term.dt <- merge(
    target.count.per.term.dt,
    all.term.sig.dt[, .(term_name, N)],
    by = "term_name"
)

target.count.per.term.dt[, `:=`(
    target_ratio = target_count / term_count,
    target_ratio_to_N_RCC4VHL = target_count / N
)]

target.count.per.term.dt[, `:=`(
    systematic_flag = target_ratio_to_N_RCC4VHL >= 0.1
)]

target.count.per.term.dt[order(target_ratio_to_N_RCC4VHL, decreasing = TRUE)]
```

```
##                            term_name      term_group             method
##  1:             Cytoplasmic ribosome Gene expression Ribosome profiling
##  2:        Pentose phosphate pathway      Metabolism Ribosome profiling
##  3:           Mitochondrial ribosome Gene expression Polysome profiling
##  4:              Translation factors Gene expression Ribosome profiling
##  5:                       Proteasome Gene expression Polysome profiling
##  6: Chaperones and folding catalysts Gene expression Polysome profiling
##  7:          Fatty acid biosynthesis      Metabolism Polysome profiling
##  8:            Amino acid metabolism      Metabolism Polysome profiling
##  9:                       Glycolysis      Metabolism Ribosome profiling
## 10:                           Oxphos      Metabolism Polysome profiling
## 11:                        TCA cycle      Metabolism Ribosome profiling
## 12:                 Ubiquitin system Gene expression Polysome profiling
## 13:         Messenger RNA biogenesis Gene expression Ribosome profiling
## 14:           Fatty acid degradation      Metabolism Polysome profiling
## 15:                      Spliceosome Gene expression Ribosome profiling
## 16:          Transcription machinery Gene expression Polysome profiling
## 17:            Nucleotide metabolism      Metabolism Ribosome profiling
## 18:            Transcription factors Gene expression Polysome profiling
## 19:             Membrane trafficking Gene expression Polysome profiling
##     target_count term_count   N target_ratio target_ratio_to_N_RCC4VHL
##  1:           71         91  75  0.780219780                0.94666667
##  2:            3         30  20  0.100000000                0.15000000
##  3:           10         77  74  0.129870130                0.13513514
##  4:           10         97  81  0.103092784                0.12345679
##  5:            5         75  60  0.066666667                0.08333333
##  6:           12        231 166  0.051948052                0.07228916
##  7:            2         44  28  0.045454545                0.07142857
##  8:           11        293 163  0.037542662                0.06748466
##  9:            2         67  39  0.029850746                0.05128205
## 10:            5        121 102  0.041322314                0.04901961
## 11:            1         30  24  0.033333333                0.04166667
## 12:           18        846 526  0.021276596                0.03422053
## 13:           11        472 329  0.023305085                0.03343465
## 14:            1         43  30  0.023255814                0.03333333
## 15:            8        411 307  0.019464720                0.02605863
## 16:            5        255 194  0.019607843                0.02577320
## 17:            2        154  84  0.012987013                0.02380952
## 18:           11       1311 486  0.008390542                0.02263374
## 19:           10       1531 898  0.006531679                0.01113586
##     systematic_flag
##  1:            TRUE
##  2:            TRUE
##  3:            TRUE
##  4:            TRUE
##  5:           FALSE
##  6:           FALSE
##  7:           FALSE
##  8:           FALSE
##  9:           FALSE
## 10:           FALSE
## 11:           FALSE
## 12:           FALSE
## 13:           FALSE
## 14:           FALSE
## 15:           FALSE
## 16:           FALSE
## 17:           FALSE
## 18:           FALSE
## 19:           FALSE
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
##  [1] Rcpp_1.0.4.6       mvtnorm_1.1-1      lattice_0.20-41    class_7.3-17      
##  [5] multcompView_0.1-8 zoo_1.8-8          digest_0.6.25      lmtest_0.9-37     
##  [9] R6_2.4.1           plyr_1.8.6         EMT_1.1            stats4_4.0.0      
## [13] evaluate_0.14      rootSolve_1.8.2.1  e1071_1.7-3        pillar_1.4.4      
## [17] rlang_0.4.10       Exact_2.1          multcomp_1.4-15    rstudioapi_0.11   
## [21] Matrix_1.2-18      labeling_0.3       splines_4.0.0      munsell_0.5.0     
## [25] compiler_4.0.0     xfun_0.14          pkgconfig_2.0.3    libcoin_1.0-6     
## [29] DescTools_0.99.38  htmltools_0.4.0    tidyselect_1.1.0   tibble_3.0.1      
## [33] lmom_2.8           expm_0.999-5       coin_1.3-1         codetools_0.2-16  
## [37] matrixStats_0.56.0 crayon_1.3.4       withr_2.4.1        rcompanion_2.3.26 
## [41] MASS_7.3-51.6      grid_4.0.0         gtable_0.3.0       lifecycle_0.2.0   
## [45] scales_1.1.1       gld_2.6.2          stringi_1.4.6      farver_2.0.3      
## [49] ellipsis_0.3.1     generics_0.0.2     vctrs_0.3.1        boot_1.3-25       
## [53] sandwich_3.0-0     nortest_1.0-4      TH.data_1.0-10     tools_4.0.0       
## [57] glue_1.4.1         purrr_0.3.4        survival_3.1-12    yaml_2.2.1        
## [61] colorspace_1.4-1   modeltools_0.2-23
```
