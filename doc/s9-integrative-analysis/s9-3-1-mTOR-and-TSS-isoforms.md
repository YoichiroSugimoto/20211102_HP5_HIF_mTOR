s9-3-1 mTOR and TSS isoforms
================
Yoichiro Sugimoto
10 February, 2022

  - [Environment setup and data
    preprocessing](#environment-setup-and-data-preprocessing)
  - [Data import](#data-import)
  - [Analysis of TSS isoforms with differential sensitivity to mTOR in
    relation to the sequence
    features](#analysis-of-tss-isoforms-with-differential-sensitivity-to-mtor-in-relation-to-the-sequence-features)
  - [Examples](#examples)
      - [Data import](#data-import-1)
      - [DERA](#dera)
      - [Other examples](#other-examples)
  - [Session information](#session-information)

# Environment setup and data preprocessing

``` r
## Specify the number of CPUs to be used
processors <- 8

## library("BiocParallel")
## register(MulticoreParam(processors))

temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)  
source(file.path("../s6-differential-expression-and-tss-usage/functions/load_total_analysis_results.R"), chdir = TRUE)
```

    ## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
    ## [1] "The following R objects were exported: total.sample.dt, total.coldata.df, total.comparison.dt"
    ## [1] "Comparison information was loaded"
    ## [1] "/camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/results"
    ## [1] "The following objects were loaded: tss.de.res.dt, tss.ratio.res.dt, diff.tss.res.dt"

``` r
source(file.path("../s8-analysis-of-translation/functions/test_differential_translation-v2.R"))
```

    ## [1] "The following functions were exported: analyzeDtg(), subsetColdata()"

``` r
s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.2.tx.assignment.dir <- file.path(s4.tss.dir, "s4-2-transcript-assignment")
s4.2.1.tss.tx.map.RCC4.dir <- file.path(s4.2.tx.assignment.dir, "s4-2-1-tss-transcript-mapping-RCC4")

s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs")
s8.1.1.dir <- file.path(s8.1.dir, "gene-level-dte")
s8.1.2.dir <- file.path(s8.1.dir, "tx-level-dte")
s8.2.dte.iso.dir <- file.path(s8.dir, "s8-2-differentially-translated-isoforms")
s8.3.dir <- file.path(s8.dir, "s8-3-validation-of-method")

s9.dir <- file.path(results.dir, "s9-integrative-analysis")

set.seed(0)
```

# Data import

``` r
mtor.tss.sig.dt <- file.path(s8.2.dte.iso.dir, "RCC4_noVHL_EIF4E2_(NA|Torin1).csv") %>% fread

all.filtered.tss.dt <- file.path(
    s8.3.dir,
    "filtered_tss_for_polysome_analysis.csv"
) %>%
    fread

mtor.tss.sig.dt <- mtor.tss.sig.dt[
]

print("The number of all significant genes (at least two isoforms with FDR < 0.1)")
```

    ## [1] "The number of all significant genes (at least two isoforms with FDR < 0.1)"

``` r
nrow(
    mtor.tss.sig.dt[
        gene_FDR < 0.1 & tx_FDR < 0.1
    ][, .N, by = gene_id][N > 1]
)
```

    ## [1] 835

``` r
rcc4.vhl.mrl.dt <- file.path(
    s8.1.2.dir,
    "RCC4_VHL_EIF4E2_yy_xx__Torin1_vs_NA-mean_ribosome_loading.csv"
) %>%
    fread

tx.meta.info.dt <- file.path(s8.3.dir, "processed-tx-meta-info.csv") %>%
    fread

sig.rcc4.vhl.mrl.dt <- merge(
    rcc4.vhl.mrl.dt[
       biotype == "protein_coding",
        .(tss_name, gene_id, gene_name, MRL_treated, MRL_base, MRL_log2fc)],
    mtor.tss.sig.dt[, .(tss_name, gene_FDR, tx_FDR)],
    all.x = TRUE,
    by = "tss_name"
) %>%
    merge(
        y = tx.meta.info.dt[
          , .(tss_name, TOP_motif_length, uORF_all, uORF_all_capped, cds_len)
        ],
        by = "tss_name"
    )

sig.rcc4.vhl.mrl.dt[
        tss_name %in% all.filtered.tss.dt[
                      RCC4_VHL_NA == TRUE & RCC4_VHL_Torin1 == TRUE, tss_name
                  ]
]
```

    ##                tss_name         gene_id  gene_name MRL_treated MRL_base
    ##    1: ENSG00000000003_3 ENSG00000000003     TSPAN6    2.710074 6.089431
    ##    2: ENSG00000000003_4 ENSG00000000003     TSPAN6    2.775607 5.361219
    ##    3: ENSG00000000419_1 ENSG00000000419       DPM1    2.351660 5.401415
    ##    4: ENSG00000000971_1 ENSG00000000971        CFH    3.788103 5.419373
    ##    5: ENSG00000001036_1 ENSG00000001036      FUCA2    3.264077 5.742278
    ##   ---                                                                  
    ## 9714: ENSG00000285976_3 ENSG00000285976 AL135905.2    3.207185 3.824970
    ## 9715: ENSG00000285976_4 ENSG00000285976 AL135905.2    3.127338 3.543440
    ## 9716: ENSG00000286019_1 ENSG00000286019  NOTCH2NLB    3.909967 3.301800
    ## 9717: ENSG00000286140_2 ENSG00000286140      DERPC    3.705092 4.390841
    ## 9718: ENSG00000286522_1 ENSG00000286522       H3C2    2.491972 3.811965
    ##       MRL_log2fc     gene_FDR       tx_FDR TOP_motif_length uORF_all
    ##    1: -1.1679751 0.0168265222 0.000000e+00      0.330708661 0.000000
    ##    2: -0.9497578 0.0168265222 0.000000e+00      0.356389215 0.000000
    ##    3: -1.1996580 1.0000000000           NA      0.409090909 0.000000
    ##    4: -0.5166503           NA           NA      0.425660377 1.000000
    ##    5: -0.8149482 0.0047858788 0.000000e+00      0.057842346 0.000000
    ##   ---                                                               
    ## 9714: -0.2541407 0.0002789074 1.000000e+00      0.213182286 4.000000
    ## 9715: -0.1802154 0.0002789074 8.140383e-05      0.109116022 3.387431
    ## 9716:  0.2439037           NA           NA      0.028409091 2.000000
    ## 9717: -0.2449877 1.0000000000           NA      0.001953125 3.000000
    ## 9718: -0.6132473           NA           NA      0.362886598 0.000000
    ##       uORF_all_capped cds_len
    ##    1:               0     738
    ##    2:               0     738
    ##    3:               0     783
    ##    4:               1    3696
    ##    5:               0    1404
    ##   ---                        
    ## 9714:              3+     522
    ## 9715:              3+     522
    ## 9716:               2     750
    ## 9717:              3+     399
    ## 9718:               0     411

``` r
sig.rcc4.vhl.mrl.dt <- sig.rcc4.vhl.mrl.dt[
    ## gene_FDR < 0.1 & tx_FDR < 0.1
]
```

# Analysis of TSS isoforms with differential sensitivity to mTOR in relation to the sequence features

``` r
sig.rcc4.vhl.mrl.dt[, `:=`(
    max_MRL_log2fc = max(MRL_log2fc),
    min_MRL_log2fc = min(MRL_log2fc)
), by = gene_id]

sig.rcc4.vhl.mrl.dt[, `:=`(
    TSS_mTOR_group = case_when(
        MRL_log2fc == max_MRL_log2fc ~ "Resistant",
        MRL_log2fc == min_MRL_log2fc ~ "Hypersensitive"
    ),
    dMRL_log2fc = max_MRL_log2fc - min_MRL_log2fc
)]

sl.sig.rcc4.vhl.mrl.dt <- sig.rcc4.vhl.mrl.dt[
    !is.na(TSS_mTOR_group) & abs(dMRL_log2fc) > 0
][order(gene_name, TSS_mTOR_group)]

te.diff.classes <- c(
    "Small (1 ~ 1.2)", "Medium (1.2 ~ 1.5)", "Large (> 1.5)"
)

sl.sig.rcc4.vhl.mrl.dt[, `:=`(
    te_diff_by_tss = case_when(
        abs(dMRL_log2fc) > log2(1.5) ~ "Large (> 1.5)",
        abs(dMRL_log2fc) > log2(1.2) ~ "Medium (1.2 ~ 1.5)",
        abs(dMRL_log2fc) > 0 ~ "Small (1 ~ 1.2)"
    ) %>% factor(levels = te.diff.classes)
)]

d.sl.sig.rcc4.vhl.mrl.dt <- dcast(
    sl.sig.rcc4.vhl.mrl.dt,
    gene_id + gene_name + te_diff_by_tss ~ TSS_mTOR_group,
    value.var = c(
        "tss_name", "MRL_log2fc", "TOP_motif_length", "uORF_all", "cds_len"
    )
)

print("All significant genes with annotation of multiple TSS available")
```

    ## [1] "All significant genes with annotation of multiple TSS available"

``` r
d.sl.sig.rcc4.vhl.mrl.dt[, table(te_diff_by_tss) %>% addmargins]
```

    ## te_diff_by_tss
    ##    Small (1 ~ 1.2) Medium (1.2 ~ 1.5)      Large (> 1.5)                Sum 
    ##               2039               1574                977               4590

``` r
## Test function
runWilcox <- function(te.diff.class, d.sl.sig.rcc4.vhl.mrl.dt, sl.genes, test.col){
    all.dt <- d.sl.sig.rcc4.vhl.mrl.dt[gene_id %in% sl.genes] 
    
    test.dt <- all.dt[te_diff_by_tss == te.diff.class] 
    
    wil.p <-  test.dt %>%
    {wilcox.test(
         .[, get(paste0(test.col, "_Hypersensitive"))],
         .[, get(paste0(test.col, "_Resistant"))],
         alternative = "two.sided",
         paired = TRUE
     )$p.value}

    wil.res.dt <- data.table(
        tested_data = test.col,
        te_diff_by_tss = te.diff.class,
        N = nrow(test.dt),
        all_N = nrow(all.dt),
        wilcox_p = wil.p
    )
    return(wil.res.dt)
}

sig.th <- 0.05

## TOP mptof length
print("TOP motif length")
```

    ## [1] "TOP motif length"

``` r
diff.top.genes <- d.sl.sig.rcc4.vhl.mrl.dt[
  TOP_motif_length_Hypersensitive != TOP_motif_length_Resistant, gene_id
]

top.test.res.dt <- lapply(
    te.diff.classes,
    runWilcox,
    d.sl.sig.rcc4.vhl.mrl.dt = d.sl.sig.rcc4.vhl.mrl.dt,
    sl.genes = diff.top.genes,
    test.col = "TOP_motif_length"
) %>%
    rbindlist

top.test.res.dt[, padj := p.adjust(wilcox_p, method = "holm")]
top.test.res.dt[, sig_mark := case_when(
                      padj < sig.th * 0.1 ~ "**",
                      padj < sig.th ~ "*",
                      TRUE ~ NA_character_
                  )
                ]
print(top.test.res.dt)
```

    ##         tested_data     te_diff_by_tss    N all_N     wilcox_p         padj
    ## 1: TOP_motif_length    Small (1 ~ 1.2) 2019  4535 3.429125e-12 3.429125e-12
    ## 2: TOP_motif_length Medium (1.2 ~ 1.5) 1562  4535 7.058822e-60 2.117647e-59
    ## 3: TOP_motif_length      Large (> 1.5)  954  4535 4.417573e-44 8.835146e-44
    ##    sig_mark
    ## 1:       **
    ## 2:       **
    ## 3:       **

``` r
ggplot(
    data = sl.sig.rcc4.vhl.mrl.dt[gene_id %in% diff.top.genes],
    aes(
        x = te_diff_by_tss,
        y = TOP_motif_length
    )
) +
    geom_boxplot(aes(fill = TSS_mTOR_group), outlier.shape = NA) +
    theme(aspect.ratio = 1.5) +
    scale_fill_bright(name = "Sensitivity of isoform to mTOR inhibition") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    coord_cartesian(ylim = c(0, 4.5)) +
    xlab("MRL log2 fold change difference between TSS isoforms") +
    ylab("TOP motif length")
```

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/TSS_dependent_sensitivity_to_mTOR-1.png)<!-- -->

``` r
## uORF length
print("uORF length")
```

    ## [1] "uORF length"

``` r
diff.uorf.genes <- d.sl.sig.rcc4.vhl.mrl.dt[
  uORF_all_Hypersensitive != uORF_all_Resistant, gene_id
]

uorf.test.res.dt <- lapply(
    te.diff.classes,
    runWilcox,
    d.sl.sig.rcc4.vhl.mrl.dt = d.sl.sig.rcc4.vhl.mrl.dt,
    sl.genes = diff.uorf.genes,
    test.col = "uORF_all"
) %>%
    rbindlist

uorf.test.res.dt[, padj := p.adjust(wilcox_p, method = "holm")]
uorf.test.res.dt[, sig_mark := case_when(
                      padj < sig.th * 0.1 ~ "**",
                      padj < sig.th ~ "*",
                      TRUE ~ NA_character_
                  )
                ]
print(uorf.test.res.dt)
```

    ##    tested_data     te_diff_by_tss   N all_N     wilcox_p         padj sig_mark
    ## 1:    uORF_all    Small (1 ~ 1.2) 887  2124 1.907965e-02 1.907965e-02        *
    ## 2:    uORF_all Medium (1.2 ~ 1.5) 724  2124 2.317370e-17 6.952109e-17       **
    ## 3:    uORF_all      Large (> 1.5) 513  2124 2.443771e-11 4.887542e-11       **

``` r
ggplot(
    data = sl.sig.rcc4.vhl.mrl.dt[gene_id %in% diff.uorf.genes],
    aes(
        x = te_diff_by_tss,
        y = uORF_all,
        fill = TSS_mTOR_group
    )
) +
    geom_boxplot(outlier.shape = NA) +
    theme(aspect.ratio = 1.5) +
    scale_fill_bright(name = "Sensitivity of isoform to mTOR inhibition") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    coord_cartesian(ylim = c(0, 7)) +
    xlab("MRL log2 fold change difference between TSS isoforms") +
    ylab("uORF number")
```

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/TSS_dependent_sensitivity_to_mTOR-2.png)<!-- -->

``` r
d.sl.sig.rcc4.vhl.mrl.dt[
    gene_id %in% mtor.tss.sig.dt[gene_FDR < 0.1, gene_id] &
    tss_name_Hypersensitive %in% mtor.tss.sig.dt[tx_FDR < 0.1, tss_name] &
    tss_name_Resistant %in% mtor.tss.sig.dt[tx_FDR < 0.1, tss_name] &
    te_diff_by_tss == "Large (> 1.5)", .(gene_name, gene_id, TOP_motif_length_Hypersensitive, TOP_motif_length_Resistant)
]
```

    ##     gene_name         gene_id TOP_motif_length_Hypersensitive
    ##  1:    RALBP1 ENSG00000017797                     0.002659574
    ##  2:  SERPINB1 ENSG00000021355                     1.883325883
    ##  3:      TYMP ENSG00000025708                     0.268429938
    ##  4:      NFYC ENSG00000066136                     2.693094629
    ##  5:    FAM50A ENSG00000071859                     0.058383234
    ##  6:     KEAP1 ENSG00000079999                     2.230683091
    ##  7:     PSMC5 ENSG00000087191                     0.537360390
    ##  8:      SUN2 ENSG00000100242                     1.290123457
    ##  9:       YY1 ENSG00000100811                     3.823673469
    ## 10:    PABPN1 ENSG00000100836                     0.027522936
    ## 11:     RAB2A ENSG00000104388                     0.215360664
    ## 12:     BCAT2 ENSG00000105552                     1.831683168
    ## 13:      CAV1 ENSG00000105974                     1.981703075
    ## 14:     HSPB1 ENSG00000106211                     0.074786461
    ## 15:     CRYAB ENSG00000109846                     0.024095298
    ## 16:       MVK ENSG00000110921                     0.555429864
    ## 17:      SOD2 ENSG00000112096                     1.572104019
    ## 18:      BRD8 ENSG00000112983                     2.060344828
    ## 19:     ABCB6 ENSG00000115657                     0.500000000
    ## 20:    TMEM59 ENSG00000116209                     0.035903840
    ## 21:     SERP1 ENSG00000120742                     4.238884045
    ## 22:      SSPN ENSG00000123096                     0.105555556
    ## 23:    WRNIP1 ENSG00000124535                     0.783001808
    ## 24:    GTF2F1 ENSG00000125651                     3.703543022
    ## 25:      ASS1 ENSG00000130707                     1.881351039
    ## 26:     MAP1B ENSG00000131711                     2.240506329
    ## 27:      RARA ENSG00000131759                     0.032258065
    ## 28:    LRRC41 ENSG00000132128                     1.544943820
    ## 29:     MUTYH ENSG00000132781                     6.646288210
    ## 30:     TPGS2 ENSG00000134779                     0.592749929
    ## 31:     CKAP4 ENSG00000136026                     2.020552344
    ## 32:      ATIC ENSG00000138363                     0.480685921
    ## 33:    SCARB2 ENSG00000138760                     3.393442623
    ## 34:     RITA1 ENSG00000139405                     0.079831933
    ## 35:     RAB20 ENSG00000139832                     1.073943662
    ## 36:    SEC11A ENSG00000140612                     0.301401869
    ## 37:     G6PC3 ENSG00000141349                     0.225261217
    ## 38:   SLC44A3 ENSG00000143036                     0.921475875
    ## 39:      RALB ENSG00000144118                     0.045170257
    ## 40:     SNX12 ENSG00000147164                     2.385146805
    ## 41:   TMSB15B ENSG00000158427                     2.240000000
    ## 42:      SHC1 ENSG00000160691                     0.031545741
    ## 43:  PDZK1IP1 ENSG00000162366                     5.273309766
    ## 44:   TM4SF18 ENSG00000163762                     0.154900617
    ## 45:    GPR146 ENSG00000164849                     0.176470588
    ## 46:      NEMF ENSG00000165525                     0.085714286
    ## 47:    CACNB3 ENSG00000167535                     6.069444444
    ## 48:    TUBA1A ENSG00000167552                     2.919431280
    ## 49:    IGFBP6 ENSG00000167779                     0.094213650
    ## 50:     MLST8 ENSG00000167965                     0.542386185
    ## 51:     CDC40 ENSG00000168438                     6.142857143
    ## 52:     PGAM1 ENSG00000171314                     0.102315161
    ## 53:     RMDN1 ENSG00000176623                     0.713469388
    ## 54:     YIPF6 ENSG00000181704                     0.022185247
    ## 55:     GPR19 ENSG00000183150                     0.900900901
    ## 56:      APOO ENSG00000184831                     0.500000000
    ## 57:    PLSCR3 ENSG00000187838                     1.343465046
    ## 58:      PJA2 ENSG00000198961                     0.330555556
    ## 59:    TRIM16 ENSG00000221926                     0.000000000
    ## 60:     THTPA ENSG00000259431                     3.204545455
    ##     gene_name         gene_id TOP_motif_length_Hypersensitive
    ##     TOP_motif_length_Resistant
    ##  1:                0.298245614
    ##  2:                0.229910714
    ##  3:                0.573369565
    ##  4:                0.063218391
    ##  5:                2.341772152
    ##  6:                0.591379310
    ##  7:                0.003021148
    ##  8:                0.044943820
    ##  9:                0.018181818
    ## 10:                0.124031008
    ## 11:                0.115432742
    ## 12:                0.348993289
    ## 13:                0.717245358
    ## 14:                0.102772643
    ## 15:                0.329938900
    ## 16:                0.369918699
    ## 17:                1.142857143
    ## 18:                0.080903104
    ## 19:                0.387755102
    ## 20:                1.219954649
    ## 21:                0.054278416
    ## 22:                0.169291339
    ## 23:                0.208333333
    ## 24:                0.050387597
    ## 25:                0.463576159
    ## 26:                0.307500889
    ## 27:                0.610486891
    ## 28:                0.033664881
    ## 29:                0.281250000
    ## 30:                0.344942571
    ## 31:                0.064285714
    ## 32:                0.271914132
    ## 33:                0.332268371
    ## 34:                0.040983607
    ## 35:                0.010000000
    ## 36:                0.161712247
    ## 37:                0.048442907
    ## 38:                0.329305136
    ## 39:                0.188191882
    ## 40:                0.407608696
    ## 41:                0.230769231
    ## 42:                0.108306189
    ## 43:                3.276269185
    ## 44:                0.122580645
    ## 45:                0.000000000
    ## 46:                0.091447368
    ## 47:                0.030805687
    ## 48:                0.012275618
    ## 49:                0.032116788
    ## 50:                0.941721854
    ## 51:                0.054708520
    ## 52:                0.028985507
    ## 53:                0.000000000
    ## 54:                0.034883721
    ## 55:                0.555555556
    ## 56:                0.181969950
    ## 57:                0.061855670
    ## 58:                0.038486628
    ## 59:                0.358565737
    ## 60:                0.150000000
    ##     TOP_motif_length_Resistant

# Examples

## Data import

``` r
r4vhl.count.per.fraction.dt <- fread(file.path(
    s8.1.2.dir,
    "RCC4_VHL_EIF4E2_yy_xx__Torin1_vs_NA-normalized_count.csv"
))

r4vhl.prop.dt <- proportionPerFraction(
    r4vhl.count.per.fraction.dt, ref.col = "tss_name", data.col.grep = "^RCC4"
)
r4vhl.prop.dt[, `:=`(
    tss_index = str_split_fixed(tss_name, "_", 2)[, 2],
    fraction2 = gsub("ribo", "", fraction) %>% {gsub("8", "8+", .)}
)]

r4novhl.count.per.fraction.dt <- fread(file.path(
    s8.1.2.dir,
    "RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA-normalized_count.csv"
))

r4novhl.prop.dt <- proportionPerFraction(
    r4novhl.count.per.fraction.dt, ref.col = "tss_name", data.col.grep = "^RCC4"
)
r4novhl.prop.dt[, `:=`(
    tss_index = str_split_fixed(tss_name, "_", 2)[, 2],
    fraction2 = gsub("ribo", "", fraction) %>% {gsub("8", "8+", .)}
)]
```

## DERA

DERA (ENSG00000023697)

``` r
treatment.colors <- c(
    "NA" = "gray40",
    "Torin1" = "darkslateblue"
)

## mTOR sensitivity of each isoform
ggplot(
    data = r4novhl.prop.dt[gene_id == "ENSG00000023697"],
    aes(
        x = fraction2,
        y = mean_ratio,
        color = treatment,
        group = treatment
    )
) +
    geom_ribbon(
        aes(
            ymin = lower_ratio_range, ymax = upper_ratio_range, fill = treatment
        ),
        color = NA, alpha = 0.2
    ) +
    geom_line(
        size = 1.25
    ) +
    ylab("% of mRNA in each polysome fraction in RCC-4") +
    xlab("Fraction") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = treatment.colors) +
    scale_fill_manual(values = treatment.colors) +
    coord_cartesian(ylim = c(0, 0.4)) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 2,
        legend.position = "bottom"
    ) +
    facet_grid(~ tss_index)
```

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/DERA-1.png)<!-- -->

``` r
## mTOR sensitivity by isoform and fraction
ggplot(
    data = rbind(
        r4novhl.prop.dt[gene_id == "ENSG00000023697"],
        r4vhl.prop.dt[gene_id == "ENSG00000023697" & treatment == "NA"]
    ),
    aes(
        x = fraction2,
        y = mean_ratio_by_fraction,
        color = tss_index,
        group = tss_index
    )
) +
    geom_ribbon(
        aes(
            ymin = lower_ratio_range_by_fraction,
            ymax = upper_ratio_range_by_fraction,
            fill = tss_index
        ),
        color = NA, alpha = 0.2
    ) +
    geom_line(
        size = 1.25
    ) +
    ylab("% of mRNA isoform by each polysome fraction") +
    xlab("Fraction") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(0, 0.9)) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 2,
        legend.position = "bottom"
    ) +
    facet_grid(~ factor(
                   paste(VHL, treatment),
                   levels = c("VHL NA", "noVHL NA", "noVHL Torin1")
               ))
```

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/DERA-2.png)<!-- -->

## Other examples

Other examples

SERP1 (ENSG00000120742); stress associated endoplasmic reticulum protein
SCARB2 (ENSG00000138760); Scavenger receptor class B member 2 GTF2F1
(ENSG00000125651); General transcription factor IIF subunit 1

``` r
## To select top3 abundanct isoforms
total.rcc4.vhl.loss.comp.name <- "RCC4_xx_HIF1B_N__noVHL_vs_VHL"
rcc4.tss.de.res.dt <- tss.de.res.dt[
    comparison_name == total.rcc4.vhl.loss.comp.name,
    .(tss_name, gene_id, meanNormCount_base)
][tss_name %in% r4vhl.prop.dt[, tss_name]]

top3.tss.dt <- rcc4.tss.de.res.dt[
    order(gene_id, meanNormCount_base, decreasing = TRUE)
][, head(.SD, 3), by = list(gene_id)]

r4vhl.prop.dt <- r4vhl.prop.dt[tss_name %in% top3.tss.dt[, tss_name]]

sl.genes <- c(
    "SERP1" = "ENSG00000120742",
    "SCARB2" = "ENSG00000138760",
    "GTF2F1" = "ENSG00000125651"
)

plotTssExamples <- function(sl.gene.id, sl.gene.name, r4vhl.prop.dt, tx.meta.info.dt){

    tx.meta.info.dt[
        tss_name %in% r4vhl.prop.dt[gene_id == sl.gene.id, tss_name],
        .(
            tss_name,
            TOP_motif_length = round(TOP_motif_length, digits = 1),
            uORF_all = round(uORF_all, digits = 1)
        )
    ] %>% print

    g1 <- ggplot(
        data = r4vhl.prop.dt[gene_id == sl.gene.id],
        aes(
            x = fraction2,
            y = mean_ratio,
            color = treatment,
            group = treatment
        )
    ) +
        geom_ribbon(
            aes(
                ymin = lower_ratio_range, ymax = upper_ratio_range, fill = treatment
            ),
            color = NA, alpha = 0.2
        ) +
        geom_line(
            size = 1.25
        ) +
        ylab("% of mRNA in each polysome fraction in RCC-4") +
        xlab("Fraction") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        geom_hline(yintercept = 0) +
        scale_color_manual(values = treatment.colors) +
        scale_fill_manual(values = treatment.colors) +
        coord_cartesian(ylim = c(0, 0.4)) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            aspect.ratio = 2,
            legend.position = "bottom"
        ) +
        facet_grid(~ tss_index) +
        ggtitle(sl.gene.name)

    print(g1)
    return()
}

temp <- mapply(
    plotTssExamples,
    sl.gene.id = sl.genes,
    sl.gene.name = names(sl.genes),
    list(r4vhl.prop.dt),
    list(tx.meta.info.dt)
)
```

    ##             tss_name TOP_motif_length uORF_all
    ## 1: ENSG00000120742_1              0.1        0
    ## 2: ENSG00000120742_2              0.6        0
    ## 3: ENSG00000120742_3              4.2        0

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/other_examples-1.png)<!-- -->

    ##             tss_name TOP_motif_length uORF_all
    ## 1: ENSG00000138760_2              0.2        0
    ## 2: ENSG00000138760_3              0.3        0
    ## 3: ENSG00000138760_5              3.4        0

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/other_examples-2.png)<!-- -->

    ##             tss_name TOP_motif_length uORF_all
    ## 1: ENSG00000125651_1              0.7      0.7
    ## 2: ENSG00000125651_2              0.1      0.0
    ## 3: ENSG00000125651_3              3.7      0.0

![](s9-3-1-mTOR-and-TSS-isoforms_files/figure-gfm/other_examples-3.png)<!-- -->

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.28        stringr_1.4.0     magrittr_1.5      data.table_1.12.8
    ## [5] dplyr_1.0.0       khroma_1.3.0      ggplot2_3.3.1     rmarkdown_2.2    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6     pillar_1.4.4     compiler_4.0.0   tools_4.0.0     
    ##  [5] digest_0.6.25    evaluate_0.14    lifecycle_0.2.0  tibble_3.0.1    
    ##  [9] gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.10     yaml_2.2.1      
    ## [13] xfun_0.14        withr_2.4.1      generics_0.0.2   vctrs_0.3.1     
    ## [17] grid_4.0.0       tidyselect_1.1.0 glue_1.4.1       R6_2.4.1        
    ## [21] purrr_0.3.4      farver_2.0.3     scales_1.1.1     ellipsis_0.3.1  
    ## [25] htmltools_0.4.0  colorspace_1.4-1 labeling_0.3     stringi_1.4.6   
    ## [29] munsell_0.5.0    crayon_1.3.4
