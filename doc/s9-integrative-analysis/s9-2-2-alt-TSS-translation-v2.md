---
title: "s9-2-2 The effect of alternate TSS usage on translation"
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


# Environment setup and data preprocessing


```r
## Specify the number of CPUs to be used
processors <- 8

## library("BiocParallel")
## register(MulticoreParam(processors))

temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)  
source(file.path("../s6-differential-expression-and-tss-usage/functions/load_total_analysis_results.R"), chdir = TRUE)
```

```
## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
## [1] "The following R objects were exported: total.sample.dt, total.coldata.df, total.comparison.dt"
## [1] "Comparison information was loaded"
## [1] "/camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/results"
## [1] "The following objects were loaded: tss.de.res.dt, tss.ratio.res.dt, diff.tss.res.dt"
```

```r
source(file.path("../s8-analysis-of-translation/functions/test_differential_translation-v2.R"))
```

```
## [1] "The following functions were exported: analyzeDtg(), subsetColdata()"
```

```r
s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs")
s8.1.1.dir <- file.path(s8.1.dir, "gene-level-dte")
s8.1.2.dir <- file.path(s8.1.dir, "tx-level-dte")
s8.2.dte.iso.dir <- file.path(s8.dir, "s8-2-differentially-translated-isoforms")
s8.3.dir <- file.path(s8.dir, "s8-3-validation-of-method")

s9.dir <- file.path(results.dir, "s9-integrative-analysis")

set.seed(0)
```

# Data preparation

## Data import

Alternate TSS usage data will be imported.


```r
## Alternative TSS usage data of all alternative TSS genes
sl.tss.all.res.dt <- fread(
    file.path(
        s9.dir,
        "alternative-vs-base-TSS-expression-regulation.csv"
    )
)

## Translation change data of RCC4 VHL loss
## Tx level
dte.res.dt <- fread(
    file.path(
        s8.1.2.dir,
        "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL.csv"
    )
)
 
## Gene level
dte.gene.res.dt <- fread(
    file.path(
        s8.1.1.dir,
        "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL.csv"
    )
)

setnames(
    dte.gene.res.dt,
    old = c("MRL_log2fc", "MRL_treated", "MRL_base"),
    new = c("gene_MRL_log2fc", "gene_MRL_treated", "gene_MRL_base")
)


## Isoform dependent translation efficiency difference test
dte.iso.dt <- fread(
    file.path(
        s8.2.dte.iso.dir,
        "RCC4_noVHL_EIF4E2_NA.csv"
    )
)

setnames(
    dte.iso.dt,
    old = c("gene_FDR", "tx_FDR"),
    new = c("gene_FDR_for_isoDTE", "tx_FDR_for_isoDTE")
)

sl.tss.all.trsl.res.dt <- merge(
    sl.tss.all.res.dt,
    dte.gene.res.dt[, .(
        gene_id,
        padj_translation, translational_regulation,
        gene_MRL_log2fc, gene_MRL_treated, gene_MRL_base
    )],
    by = "gene_id"
) %>%
    merge(
        y = dte.res.dt[
            , .(tss_name, MRL_log2fc, MRL_treated, MRL_base)
        ],
        by = "tss_name"
    ) %>%
    merge(
        y = dte.iso.dt[!duplicated(gene_id), .(gene_id, gene_FDR_for_isoDTE)],
        all.x = TRUE,
        by = "gene_id"
    ) %>%
    merge(
        y = dte.iso.dt[, .(tss_name, tx_FDR_for_isoDTE)],
        all.x = TRUE,
        by = "tss_name"
    )


sl.tss.all.trsl.res.dt[, `:=`(
    proportion_treated_sum = sum(proportion_treated),
    proportion_base_sum = sum(proportion_base)
), by = gene_id]

sl.tss.all.trsl.res.dt[, `:=`(
    corrected_proportion_treated = proportion_treated / proportion_treated_sum,
    corrected_proportion_base = proportion_base / proportion_base_sum
)]

sl.tss.all.trsl.res.dt[, `:=`(
    dProportion = corrected_proportion_treated - corrected_proportion_base
)]

fwrite(
    sl.tss.all.trsl.res.dt,
    file.path(
        s9.dir,
        "VHL-dependent-alternate-TSS-and-translation-long.csv"
    )    
)
```


### Data for filtration



```r
all.filtered.tss.dt <- fread(
    file.path(
        s8.3.dir,
        "filtered_tss_for_polysome_analysis.csv"
    )
)
```


# Translational efficiency differences between alternate TSS isoforms against other isoforms of same gene



```r
alt.tss.res.dt <- sl.tss.all.trsl.res.dt[base_alt_TSS_flag == "alternative_TSS"]
alt.tss.res.dt[, `:=`(
    alt_TSS_diff_poly = if_else(
        gene_FDR_for_isoDTE < 0.1 & tx_FDR_for_isoDTE < 0.1, "Different", "Not different",
        missing = "Not different"
    )
)]

print(paste0(
    "The number of alternative TSS genes before filteration: ",
    nrow(alt.tss.res.dt)
))
```

```
## [1] "The number of alternative TSS genes before filteration: 149"
```

```r
print("The number of genes after filtration for translational efficiency analysis:")
```

```
## [1] "The number of genes after filtration for translational efficiency analysis:"
```

```r
alt.tss.res.dt <- alt.tss.res.dt[
    tss_name %in% all.filtered.tss.dt[RCC4_noVHL_NA == TRUE, tss_name]
]

print("The number and proportion of genes where the alternate TSS mRNA isoforms had different polysome distribution compared to other isoforms of same gene")
```

```
## [1] "The number and proportion of genes where the alternate TSS mRNA isoforms had different polysome distribution compared to other isoforms of same gene"
```

```r
alt.tss.res.dt[, table(alt_TSS_diff_poly)] %T>%
    print %>%
    prop.table %>%
    round(digit = 2)
```

```
## alt_TSS_diff_poly
##     Different Not different 
##            75            42
```

```
## alt_TSS_diff_poly
##     Different Not different 
##          0.64          0.36
```



```r
## From here MRL of isoforms from alternate TSS and other isoforms will be calculated
normalized.polysome.count.dt <- file.path(
    s8.1.2.dir, "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL-normalized_count.csv"
) %>%
    fread

ref.column.name <- "tss_name"
base.input.names <- c("RCC4_noVHL", "RCC4_VHL")

sample.colnames <- grep("^RCC4", colnames(normalized.polysome.count.dt), value = TRUE)
normalized.polysome.count.dt <- normalized.polysome.count.dt[
  , c(ref.column.name, sample.colnames), with = FALSE
]

## Classify alternate vs other TSS
normalized.polysome.count.dt <- normalized.polysome.count.dt[
    tss_name %in% sl.tss.all.res.dt[, tss_name]
]

normalized.polysome.count.dt[, `:=`(
    gene_id = str_split_fixed(tss_name, "_", 2)[, 1],
    alt_tss_flag =
        tss_name %in% sl.tss.all.res.dt[
                          base_alt_TSS_flag == "alternative_TSS", tss_name
                      ]
)]

normalized.polysome.count.dt[, `:=`(
    tss_name = if_else(
        alt_tss_flag == TRUE, paste0(gene_id, "_alternate"), paste0(gene_id, "_others")
    ),
    gene_id = NULL,
    alt_tss_flag = NULL
)]

cl.normalized.polysome.count.dt <- normalized.polysome.count.dt[
    , lapply(.SD, sum), by = tss_name,
    .SDcols = sample.colnames
]

sl.ribo.sample.groups <- grep(
    "ribo1", sample.colnames, value = TRUE
) %>%
    gsub("ribo1", "ribo[[:digit:]](|A|B)", .)

mrl.count.dt <- calculateWeightedMRL(
    normalized.polysome.count.dt = cl.normalized.polysome.count.dt,
    sl.ribo.sample.groups = sl.ribo.sample.groups,
    ref.column.name = ref.column.name
)

mrl.count4export.dt <- calcMrlStat(
    mrl.count.dt = mrl.count.dt,
    annot.df = mrl.count.dt[, "tss_name", with = FALSE],
    ref.column.name = ref.column.name,
    base.input.names = base.input.names
)

mrl.count4export.dt[, `:=`(
    gene_id = str_split_fixed(tss_name, "_", n = 2)[, 1],
    tss_group = str_split_fixed(tss_name, "_", n = 2)[, 2]
)]

bin.mrl.alt.other.dt <- dcast(
    mrl.count4export.dt,
    gene_id ~ tss_group,
    value.var = "MRL_treated"
)

bin.mrl.alt.other.dt[, MRL_log2fc_alt_others := log2(alternate / others)]

bin.mrl.alt.other.dt <- merge(
    alt.tss.res.dt[, .(
        gene_id, gene_name, alt_TSS_diff_poly,
        dProportion, alt_tss_reg_mode
    )],
    bin.mrl.alt.other.dt,
    by = "gene_id"
) %>% {
    .[order(
        alt_tss_reg_mode,
        if_else(
            alt_tss_reg_mode == "Up",
            MRL_log2fc_alt_others, -MRL_log2fc_alt_others
        )
    )][, `:=`(
         gene_name = factor(gene_name, levels = gene_name),
         alt_tss_reg_mode = factor(alt_tss_reg_mode, levels = c("Up", "Down"))
     )]
}
ggplot(
    data = bin.mrl.alt.other.dt,
    aes(
        y = gene_name,
        x = MRL_log2fc_alt_others,
        color = alt_TSS_diff_poly
    )
) +
    geom_vline(xintercept = 0, color = "gray60") +
    geom_point(
        aes(
            size = abs(dProportion)
        )
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.04))) +
    scale_color_manual(values = c(
                           "Different" = "black",
                           "Not different" = "gray60"
                       )) +
    facet_grid(alt_tss_reg_mode ~ ., scales = "free_y", space = "free_y") +
    theme_bw(12) +
    theme(
        legend.box = "vertical"
    ) +
    guides(
        size = guide_legend(
            title = "d% of alternate TSS usage\nupon VHL loss"
        ),
        color = guide_legend(
            title = "Diff. polysome distribution"
        )
    ) +
    ylab("Gene name") +
    xlab("MRL log2 fold change (alternate TSS mRNA isoform / other TSS mRNA isoform")
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/plot te differences-1.png)<!-- -->



# Simulate MRL change with omitting a parameter

Only genes with alternate TSS mRNA isoforms whose polysome distribution was significantly different from that of other isoforms from the same genes (FDR < 0.1) were analysed.

Two models will be considered here:

1. Omitting MRL change in each mRNA isoform
2. Omitting mRNA abundance change in each mRNA isoform



```r
sl.tss.all.trsl.res.dt[, `:=`(
    ## Upon VHL loss
    mean_MRL = rowMeans(cbind(
        MRL_treated, MRL_base
    ), na.rm = TRUE),
     mean_proportion = rowMeans(cbind(
        corrected_proportion_treated, corrected_proportion_base
    ), na.rm = TRUE),
    alt_TSS_diff_trsl =
        gene_id %in%
        alt.tss.res.dt[alt_TSS_diff_poly == "Different", gene_id]
)]


## Calculate necessary values
sl.tss.all.trsl.res.dt[, `:=`(
    ## Upon VHL loss
    ## Model 1
    dMRL_alt_TSS_dep_noVHL = mean_MRL * proportion_treated,
    dMRL_alt_TSS_dep_VHL = mean_MRL * proportion_base,
    ## Model 2
    dMRL_alt_TSS_indep_noVHL = MRL_treated * mean_proportion,
    dMRL_alt_TSS_indep_VHL = MRL_base * mean_proportion
)]
```


## Data export for publication



```r
for.export.dt <- copy(sl.tss.all.trsl.res.dt)[, .(
                         tss_name,
                         alt_tss_reg_mode, 
                         base_alt_TSS_flag,
                         de_log2fc,
                         tx_FDR,                         
                         corrected_proportion_base, corrected_proportion_treated,
                         dProportion,
                         meanNormCount_base, meanNormCount_treated,
                         alt_TSS_diff_trsl,
                         gene_FDR_for_isoDTE, tx_FDR_for_isoDTE,
                         gene_MRL_log2fc, MRL_log2fc, MRL_base, MRL_treated
                     )]

for.export.dt <- merge(
    filtered.tss.with.quantile.dt[, .(
        tss_name, gene_id, gene_name,
        chr, strand, start, end 
    )],
    for.export.dt,
    by = "tss_name"
)

for.export.dt[, `:=`(
    base_alt_TSS_flag = case_when(
        base_alt_TSS_flag == "alternative_TSS" ~ "alternate_TSS",
        TRUE ~ base_alt_TSS_flag
    )
)]

for.export.dt <- for.export.dt[order(
    - alt_tss_reg_mode,
    - alt_TSS_diff_trsl,
    - gene_MRL_log2fc
)]

for.export.dt[
  , alt_TSS_diff_trsl := case_when(
        gene_id %in% alt.tss.res.dt[, gene_id] ~ alt_TSS_diff_trsl
        ## Genes with too little alternate TSS expression for robust MRL calculation will be ignored
    )
]

for.export.dt <- for.export.dt[, .(
    gene_id, gene_name,
    tss_name,
    chr, strand, start, end,
    base_alt_TSS_flag,
    de_log2fc,
    tx_FDR,                         
    corrected_proportion_base, corrected_proportion_treated,
    dProportion,
    meanNormCount_base, meanNormCount_treated,
    gene_MRL_log2fc, MRL_log2fc, MRL_base, MRL_treated,
    alt_TSS_diff_trsl,
    gene_FDR_for_isoDTE, tx_FDR_for_isoDTE
)]

setnames(
    for.export.dt,
    old = c(
        "chr", "base_alt_TSS_flag",
        "de_log2fc",
        "tx_FDR",
        "corrected_proportion_base", "corrected_proportion_treated", "dProportion",
        "meanNormCount_base", "meanNormCount_treated",
        "alt_TSS_diff_trsl",
        "gene_FDR_for_isoDTE", "tx_FDR_for_isoDTE",
        "gene_MRL_log2fc", "MRL_log2fc",
        "MRL_base", "MRL_treated"
    ),
    new = c(
        "chromosome", "TSS class",
        "mRNA log2FC",
        "FDR (TSS usage change)",
        "TSS usage (RCC4 VHL)", "TSS usage (RCC4)", "delta TSS usage",
        "TPM (RCC4 VHL)", "TPM (RCC4)",
        "class (isoform dependent polysome distribution diff.)",
        "FDR (gene level)",
        "FDR (isoform level)",
        "MRL log2FC (gene level)", "MRL log2FC (isoform level)",
        "MRL (RCC4 VHL)", "MRL (RCC4)"
    )
)

fwrite(
    for.export.dt,
    file.path(
        s9.dir, "alteranate_TSS_and_translation_in_RCC4.csv"
    )
)
```


## Simulation



```r
## The analysis will be performed for genes with an alternate TSS isoform differentially translated compared to other isoforms of the same gene
sl.tss.all.trsl.res.dt <- sl.tss.all.trsl.res.dt[alt_TSS_diff_trsl == TRUE]

print("Confirmation")
```

```
## [1] "Confirmation"
```

```r
sl.tss.all.trsl.res.dt[, .(
    tss_name, gene_name,
    gene_MRL_log2fc, mean_MRL, mean_proportion
)] %>%
    {.[!complete.cases(.)]}
```

```
## Empty data.table (0 rows and 5 cols): tss_name,gene_name,gene_MRL_log2fc,mean_MRL,mean_proportion
```

```r
measured.vs.estimated.dt <- sl.tss.all.trsl.res.dt[
    ,
    list(
        gene_name,
        gene_MRL_log2fc,        
        ## For VHL loss
        simulated_gene_MRL_alt_TSS_dep_noVHL = sum(dMRL_alt_TSS_dep_noVHL),
        simulated_gene_MRL_alt_TSS_dep_VHL = sum(dMRL_alt_TSS_dep_VHL),
        simulated_gene_MRL_alt_TSS_indep_noVHL = sum(dMRL_alt_TSS_indep_noVHL),
        simulated_gene_MRL_alt_TSS_indep_VHL = sum(dMRL_alt_TSS_indep_VHL)
    ),
    by = gene_id
][!duplicated(gene_id)]

measured.vs.estimated.dt[, `:=`(
    simulated_gene_MRL_log2fc_alt_TSS_dep = log2(
        simulated_gene_MRL_alt_TSS_dep_noVHL / simulated_gene_MRL_alt_TSS_dep_VHL
    ),
    simulated_gene_MRL_log2fc_alt_TSS_indep = log2(
        simulated_gene_MRL_alt_TSS_indep_noVHL / simulated_gene_MRL_alt_TSS_indep_VHL
    )
)]

vhl.loss.lim.range <- c(-0.35, 0.5)

g1 <- ggplot(
    data = measured.vs.estimated.dt,
    aes(
        x = simulated_gene_MRL_log2fc_alt_TSS_dep,
        y = gene_MRL_log2fc,
    )
) +
    geom_vline(xintercept = 0, color = "gray60") +
    geom_hline(yintercept = 0, color = "gray60") +
    geom_smooth(method = "lm") +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(
                 label = ifelse(gene_name == "MXI1", gene_name, NA)
             ), direction = "y", nudge_y = 0.05, min.segment.length = 10) +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = vhl.loss.lim.range, ylim = vhl.loss.lim.range) +
    xlab("Simulated gene level MRL log2FC\n(altenative TSS depedent)") +
    ylab("Gene level MRL log2FC") +
    ggtitle("Model i")

g2 <- ggplot(
    data = measured.vs.estimated.dt,
    aes(
        x = simulated_gene_MRL_log2fc_alt_TSS_indep,
        y = gene_MRL_log2fc,
    )
) +
    geom_vline(xintercept = 0, color = "gray60") +
    geom_hline(yintercept = 0, color = "gray60") +
    geom_smooth(method = "lm") +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(
                 label = ifelse(gene_name == "MXI1", gene_name, NA)
             )) +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = vhl.loss.lim.range, ylim = vhl.loss.lim.range) +
    xlab("Estimated gene level MRL log2FC\n(alternative TSS independent)") +
    ylab("Gene level MRL log2FC") +
    ggtitle("Model ii")

cowplot::plot_grid(g1, g2)
```

```
## `geom_smooth()` using formula 'y ~ x'
```

```
## Warning: Removed 74 rows containing missing values (geom_text_repel).
```

```
## `geom_smooth()` using formula 'y ~ x'
```

```
## Warning: Removed 5 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

```
## Warning: Removed 74 rows containing missing values (geom_text_repel).
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/calculate blocked gene level MRL estimate-1.png)<!-- -->

```r
print("The number of genes analysed for model i")
```

```
## [1] "The number of genes analysed for model i"
```

```r
measured.vs.estimated.dt[
  , .(simulated_gene_MRL_log2fc_alt_TSS_dep, gene_MRL_log2fc)
] %>% {nrow(.[complete.cases(.)])}
```

```
## [1] 75
```

```r
measured.vs.estimated.dt %$%
    cor.test(
        x = simulated_gene_MRL_log2fc_alt_TSS_dep,
        y = gene_MRL_log2fc,
        method = "pearson",
        alternative = "two.sided"
    )
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  simulated_gene_MRL_log2fc_alt_TSS_dep and gene_MRL_log2fc
## t = 12.688, df = 73, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7422396 0.8890586
## sample estimates:
##       cor 
## 0.8294691
```

```r
print("The number of genes analysed for model ii")
```

```
## [1] "The number of genes analysed for model ii"
```

```r
measured.vs.estimated.dt[
  , .(simulated_gene_MRL_log2fc_alt_TSS_indep, gene_MRL_log2fc)
] %>% {nrow(.[complete.cases(.)])}
```

```
## [1] 70
```

```r
measured.vs.estimated.dt %$%
    cor.test(
        x = simulated_gene_MRL_log2fc_alt_TSS_indep,
        y = gene_MRL_log2fc,
        method = "pearson",
        alternative = "two.sided"
    )
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  simulated_gene_MRL_log2fc_alt_TSS_indep and gene_MRL_log2fc
## t = 5.3558, df = 68, p-value = 1.09e-06
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.3551665 0.6911944
## sample estimates:
##      cor 
## 0.544684
```

```r
cocor::cocor(
           formula = ~ gene_MRL_log2fc + simulated_gene_MRL_log2fc_alt_TSS_indep |
               gene_MRL_log2fc + simulated_gene_MRL_log2fc_alt_TSS_dep,
           data = measured.vs.estimated.dt,
           alternative = "two.sided",
           test = "williams1959"
       )
```

```
## 
##   Results of a comparison of two overlapping correlations based on dependent groups
## 
## Comparison between r.jk (gene_MRL_log2fc, simulated_gene_MRL_log2fc_alt_TSS_indep) = 0.5447 and r.jh (gene_MRL_log2fc, simulated_gene_MRL_log2fc_alt_TSS_dep) = 0.8
## Difference: r.jk - r.jh = -0.2553
## Related correlation: r.kh = 0.0729
## Data: measured.vs.estimated.dt: j = gene_MRL_log2fc, k = simulated_gene_MRL_log2fc_alt_TSS_indep, h = simulated_gene_MRL_log2fc_alt_TSS_dep
## Group size: n = 70
## Null hypothesis: r.jk is equal to r.jh
## Alternative hypothesis: r.jk is not equal to r.jh (two-sided)
## Alpha: 0.05
## 
## williams1959: Williams' t (1959)
##   t = -2.8114, df = 67, p-value = 0.0065
##   Null hypothesis rejected
```


## Plot range of gene level MRL log2 FC


```r
all.filtered.gene.dt <- fread(
    file.path(
        s8.3.dir,
        "filtered_gene_for_polysome_analysis.csv"
    )
)

dte.gene.res.dt[
    gene_id %in% all.filtered.gene.dt[
                      RCC4_noVHL_NA == TRUE & RCC4_VHL_NA == TRUE,
                      gene_id
                  ]
] %>%
    ggplot(
        aes(
            x = gene_MRL_log2fc,
            y = stat(count / sum(count))
        )
    ) +
    geom_vline(xintercept = 0, color = "gray60") +
    geom_histogram(
        binwidth = 0.025
    ) +
    coord_flip(xlim = vhl.loss.lim.range) +
    theme(
        aspect.ratio = 4
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab("Proportion") +
    xlab("MRL log2 fold change upon VHL loss")
```

```
## Warning: Removed 4 rows containing non-finite values (stat_bin).
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/plot a range of gene level MRL log2 fold change-1.png)<!-- -->

# MXI1


MXI1 was one of the most significantly differentially translated mRNAs upon VHL loss, and the analysis above indicated the strong contribution of the alternate TSS usage.
Here, the mode of action of the translational efficiency change will be examined.

## Gene level



```r
count.per.fraction.dt <- fread(file.path(
    s8.1.1.dir,
    "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL-normalized_count.csv"
))

base.cols <- c("gene_id", "gene_name", "biotype")

m.count.per.fraction.dt <- melt(
    count.per.fraction.dt,
    id.vars = base.cols,
    value.name = "normalized_count",
    variable.name = "sample_name"
)

m.count.per.fraction.dt[, `:=`(
    sample_group = sub("(.*?_.*?_.*?)_.*", "\\1", sample_name),
    VHL = str_split_fixed(sample_name, "_", n = 8)[, 2],
    EIF4E2 = str_split_fixed(sample_name, "_", n = 8)[, 3],
    clone = str_split_fixed(sample_name, "_", n = 8)[, 5],
    fraction = str_split_fixed(sample_name, "_", n = 8)[, 7]
)]

m.count.per.fraction.dt[, `:=`(
    sum_across_fraction = sum(normalized_count, na.rm = TRUE)
), by = list(gene_id, gene_name, sample_group, clone)]

m.count.per.fraction.dt[, `:=`(
    ratio_across_fraction = normalized_count / sum_across_fraction
)]

m.count.per.fraction.dt[, `:=`(
    sd_norm_count = sd(normalized_count),
    mean_norm_count = mean(normalized_count),
    sd_ratio = sd(ratio_across_fraction),
    mean_ratio = mean(ratio_across_fraction),
    lower_ratio_range = mean(ratio_across_fraction) - sd(ratio_across_fraction),
    upper_ratio_range = mean(ratio_across_fraction) + sd(ratio_across_fraction)
), by = list(gene_id, gene_name, sample_group, fraction)]

m.count.per.fraction.dt <- m.count.per.fraction.dt[clone == 1]

m.count.per.fraction.dt[gene_name == "MXI1"] %>%
    ggplot(
        aes(
            x = fraction,
            y = mean_ratio,
            color = VHL,
            group = VHL
        )
    ) +
    geom_ribbon(
        aes(
            ymin = lower_ratio_range, ymax = upper_ratio_range, fill = VHL
        ),
        color = NA, alpha = 0.2
    ) +
    geom_line(
        size = 1.25
    ) +
    ylab("% of mRNA in each polysome fraction in RCC-4") +
    xlab("Fraction") +
    ggsci::scale_color_npg() +
    ggsci::scale_fill_npg() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(0, 0.4)) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 2
    )
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/gene level-1.png)<!-- -->


## mRNA TSS isoforms


```r
count.per.fraction.dt <- fread(file.path(
    s8.1.2.dir,
    "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL-normalized_count.csv"
))

base.cols <- c("tss_name", "gene_id", "gene_name", "biotype")

m.count.per.fraction.dt <- melt(
    count.per.fraction.dt,
    id.vars = base.cols,
    value.name = "normalized_count",
    variable.name = "sample_name"
)

## Filter TSS not expressed
m.count.per.fraction.dt <- m.count.per.fraction.dt[
    tss_name %in% sl.tss.all.trsl.res.dt[, tss_name]
]

top3.tss.dt <- sl.tss.all.trsl.res.dt[
  , rowMeans(.SD), .SDcols = c("proportion_treated", "proportion_base"),
    by = list(gene_id, gene_name, tss_name)
][order(V1, decreasing = TRUE)][, head(.SD, 3), by = list(gene_id)]

m.count.per.fraction.dt[
  , tss_name := case_when(
        tss_name %in% top3.tss.dt[, tss_name] ~ tss_name,
        TRUE ~ paste0(gene_id, "_", "others")
    )
]

## Only plot top 3
m.count.per.fraction.dt <- m.count.per.fraction.dt[!grepl("_others$", tss_name)]

m.count.per.fraction.dt <- m.count.per.fraction.dt[, list(
    gene_id, gene_name, biotype, normalized_count = sum(normalized_count)
), by = list(sample_name, tss_name)][!duplicated(paste0(sample_name, "_", tss_name))]

m.count.per.fraction.dt[, `:=`(
    sample_group = sub("(.*?_.*?_.*?)_.*", "\\1", sample_name),
    VHL = str_split_fixed(sample_name, "_", n = 8)[, 2],
    EIF4E2 = str_split_fixed(sample_name, "_", n = 8)[, 3],
    clone = str_split_fixed(sample_name, "_", n = 8)[, 5],
    fraction = str_split_fixed(sample_name, "_", n = 8)[, 7],
    tss_index = str_split_fixed(tss_name, "_", n = 2)[, 2]
)]

m.count.per.fraction.dt[, `:=`(
    sum_across_fraction = sum(normalized_count, na.rm = TRUE)
), by = list(tss_name, sample_group, clone)]

m.count.per.fraction.dt[, `:=`(
    ratio_across_fraction = normalized_count / sum_across_fraction
)]

m.count.per.fraction.dt[, `:=`(
    sd_norm_count = sd(normalized_count),
    mean_norm_count = mean(normalized_count),
    sd_ratio = sd(ratio_across_fraction),
    mean_ratio = mean(ratio_across_fraction),
    lower_ratio_range = mean(ratio_across_fraction) - sd(ratio_across_fraction),
    upper_ratio_range = mean(ratio_across_fraction) + sd(ratio_across_fraction)
), by = list(tss_name, sample_group, fraction)]

m.count.per.fraction.dt <- m.count.per.fraction.dt[clone == 1]

## DE data
total.rcc4.vhl.loss.comp.name <- "RCC4_xx_HIF1B_N__noVHL_vs_VHL"
rcc4.tss.de.res.dt <- tss.de.res.dt[comparison_name == total.rcc4.vhl.loss.comp.name]

m.rcc4.tss.de.res.dt <- melt(
    rcc4.tss.de.res.dt,
    id.vars = c("tss_name", "gene_id", "gene_name"),
    measure.vars = c("meanNormCount_treated", "meanNormCount_base"),
    value.name = "TPM",
    variable.name = "data_source"
)

m.rcc4.tss.de.res.dt[
  , tss_name := case_when(
        tss_name %in% top3.tss.dt[, tss_name] ~ tss_name,
        TRUE ~ paste0(gene_id, "_", "others")
    )
]

## Only plot top 3
m.rcc4.tss.de.res.dt <- m.rcc4.tss.de.res.dt[!grepl("_others$", tss_name)]

m.rcc4.tss.de.res.dt <- m.rcc4.tss.de.res.dt[, list(
    gene_id, gene_name, TPM = sum(TPM)
), by = list(tss_name, data_source)][!duplicated(paste0(data_source, "_", tss_name))]


m.rcc4.tss.de.res.dt[, `:=`(
    VHL = case_when(
        data_source == "meanNormCount_treated" ~ "noVHL",
        data_source == "meanNormCount_base" ~ "VHL"
    ) %>% factor(levels = c("VHL", "noVHL")),
    tss_index = str_split_fixed(tss_name, "_", n = 2)[, 2]
)]


sl.gene.name <- "MXI1"

tss.to.plot.dt <- data.table(
    tss_name = unique(m.rcc4.tss.de.res.dt[gene_name == sl.gene.name, tss_name])
)[order(tss_name)]

tss.to.plot.dt[, `:=`(
    tss_index = str_split_fixed(tss_name, "_", n = 2)[, 2],
    assigned_color = c("#D33F6A", "#023FA5", "#7D87B9")#, "gray40")
)]

tss.exp.dt <- m.rcc4.tss.de.res.dt[gene_name == sl.gene.name]

tss.poly.dt <- m.count.per.fraction.dt[
    gene_name == sl.gene.name &
    VHL == "noVHL"    
]

ggplot(
    data = tss.exp.dt,
    aes(
        x = VHL,
        y = TPM,
        fill = tss_index
    )
) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
        values = setNames(tss.to.plot.dt[, assigned_color], nm = tss.to.plot.dt[, tss_index])
    ) +
    xlab("VHL status") +
    ylab("mRNA abundance (TPM)") +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 4
    )
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/MXI1 isoforms-1.png)<!-- -->

```r
ggplot(
    data = tss.poly.dt,
    aes(
        x = fraction,
        y = mean_ratio,
        color = tss_index,
        group = tss_index
    )
) +
    geom_ribbon(
        aes(
            ymin = lower_ratio_range, ymax = upper_ratio_range, fill = tss_index
        ),
        color = NA, alpha = 0.2
    ) +
    geom_line(
        size = 1.25
    ) +
    ylab("% of mRNA in each polysome fraction in RCC-4") +
    xlab("Fraction") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(
        values = setNames(tss.to.plot.dt[, assigned_color], nm = tss.to.plot.dt[, tss_index])
    ) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(0, 0.4)) +
    scale_color_manual(
        values = setNames(tss.to.plot.dt[, assigned_color], nm = tss.to.plot.dt[, tss_index])
    ) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 2
    )
```

![](s9-2-2-alt-TSS-translation-v2_files/figure-html/MXI1 isoforms-2.png)<!-- -->




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
##  [1] Rcpp_1.0.4.6     pillar_1.4.4     compiler_4.0.0   tools_4.0.0     
##  [5] digest_0.6.25    lattice_0.20-41  nlme_3.1-148     evaluate_0.14   
##  [9] lifecycle_0.2.0  tibble_3.0.1     gtable_0.3.0     mgcv_1.8-31     
## [13] pkgconfig_2.0.3  rlang_0.4.10     Matrix_1.2-18    ggsci_2.9       
## [17] ggrepel_0.8.2    yaml_2.2.1       xfun_0.14        cocor_1.1-3     
## [21] withr_2.4.1      generics_0.0.2   vctrs_0.3.1      cowplot_1.0.0   
## [25] grid_4.0.0       tidyselect_1.1.0 glue_1.4.1       R6_2.4.1        
## [29] purrr_0.3.4      farver_2.0.3     splines_4.0.0    scales_1.1.1    
## [33] ellipsis_0.3.1   htmltools_0.4.0  colorspace_1.4-1 labeling_0.3    
## [37] stringi_1.4.6    munsell_0.5.0    crayon_1.3.4
```