s8-2 Identification of differentially translated mRNA isoforms
================
Yoichiro Sugimoto
02 March, 2022

  - [Overview](#overview)
  - [Read raw and normalized count
    table](#read-raw-and-normalized-count-table)
  - [Data pre-processing](#data-pre-processing)
  - [Analysis of differentially translated mRNA
    isoforms](#analysis-of-differentially-translated-mrna-isoforms)
      - [data setup for DEXSeq
        analysis](#data-setup-for-dexseq-analysis)
      - [Define DTE analysis function](#define-dte-analysis-function)
  - [Identification of differentially translated mRNA isoforms within a
    gene](#identification-of-differentially-translated-mrna-isoforms-within-a-gene)
      - [RCC4 VHL](#rcc4-vhl)
      - [RCC4 noVHL](#rcc4-novhl)
  - [Session information](#session-information)

# Overview

Significantly differentially translated mRNA isoforms compared to other
isoforms of same gene will be identified.

``` r
## Bioconductor
library("DRIMSeq")
library("stageR")
library("DEXSeq")
## Specify the number of CPUs to be used
processors <- 24
library("BiocParallel")
BPPARAM = MulticoreParam(processors)

temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source, chdir = TRUE)
```

    ## [1] "Sample file used: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/20211102_HP5_HIF_mTOR/data/sample_data/processed_sample_file.csv"
    ## [1] "The following objects are exported: poly.coldata.df, poly.sample.dt, translation.comparison.dt"
    ## [1] "In translation.comparison.dt, xx specifies the factor compared where the comparison is specified after __, while yy is a wildcard. From left, each factor specifies cell, VHL, EIF4E2, clone, and treatment"
    ## [1] "The following functions were exported: analyzeDtg(), subsetColdata()"

``` r
set.seed(0)
```

``` r
sample.file <- file.path("../../data/sample_data/20190812_PM18011_1_sample-data.csv")
annot.dir <- normalizePath(file.path("../../annotation/"))
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)

results.dir <- file.path("../../results")

s3.alignment.stats.dir <- file.path(results.dir, "s3-alignment-statistics")
s3.4.poly.size.factor.dir <- file.path(s3.alignment.stats.dir, "polysome_size_factor")

s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.7.count.per.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-7-count-per-tss")

s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
s8.2.dte.iso.dir <- file.path(s8.dir, "s8-2-differentially-translated-isoforms")

create.dirs(c(
    s8.2.dte.iso.dir
))
```

# Read raw and normalized count table

``` r
count.per.tss.file <- file.path(s4.1.7.count.per.tss.dir, "count-per-confident-tss.csv")
count.per.tss.dt <- fread(count.per.tss.file)

## Only polysome proflied sequence data are eported
poly.count.per.tss.dt <- count.per.tss.dt[
  , c("tss_name", grep("polysome_", colnames(count.per.tss.dt), value = TRUE)),
    with = FALSE
]
colnames(poly.count.per.tss.dt) <- colnames(poly.count.per.tss.dt) %>%
    gsub("polysome_", "", .)
poly.count.per.tss.dt[, gene_id := str_split_fixed(tss_name, "_", n = 2)[, 1]]

## Only analyze transcripts from genes in primary.tx.dt
poly.count.per.tss.dt <- poly.count.per.tss.dt[
    gene_id %in% unique(primary.tx.dt[, gene_id])
]

## Read normalization factor
poly.sizefactor.dt <- fread(
    file = file.path(
        s3.4.poly.size.factor.dir,
        "library_size_factor_by_ERCC.csv"
    )
)
```

# Data pre-processing

``` r
annot.cols <- c("gene_id", "gene_name", "biotype")

poly.count.per.tss.with.meta.dt <- merge(
  primary.tx.dt[!duplicated(gene_id), annot.cols, with = FALSE],
  poly.count.per.tss.dt[, c("tss_name", "gene_id", rownames(poly.coldata.df)), with = FALSE],
  by = "gene_id",
  all.y = TRUE
)

## Ignore free fraction
poly.count.per.tss.with.meta.dt <- poly.count.per.tss.with.meta.dt[
    , !grepl("ribo0", colnames(poly.count.per.tss.with.meta.dt)), with = FALSE
]
```

# Analysis of differentially translated mRNA isoforms

## data setup for DEXSeq analysis

``` r
countDt2Df <- function(count.dt, annot.cols){
    count.df <- as.data.frame(count.dt)
    rownames(count.df) <- count.df[, "tss_name"]
    annot.df <- count.df[, c("tss_name", annot.cols)]
    count.df <- dropColumnDf(count.df, drop.vec = c("tss_name", annot.cols))
    return(list(count.df = count.df, annot.df = annot.df))
}

deseq2.poly.in.list <- countDt2Df(poly.count.per.tss.with.meta.dt, annot.cols = annot.cols)

poly.coldata.df <- poly.coldata.df[colnames(deseq2.poly.in.list$count.df), ]

if(!all(colnames(deseq2.poly.in.list$count.df) == rownames(poly.coldata.df))){
    stop("colnames of count.df does not match with rownames of poly.coldata.df")
} else {"OK"}
```

    ## [1] "OK"

``` r
annot.dt <- data.table(deseq2.poly.in.list$annot.df)
setkey(annot.dt, tss_name)
```

## Define DTE analysis function

``` r
analyzeDtIsoforms <- function(sample.base.name,
                              count.dt,
                              annot.dt,
                              input.sample.data.df,
                              dexseq.exp.design,
                              s8.2.dte.iso.dir,
                              ref.column.name = "tss_name",
                              BPPARAM){

    input.sample.data.df <-
        input.sample.data.df[
            grepl("^ribo[[:digit:]]$", input.sample.data.df[, "fraction"]),
            ]
    
    print("List of samples analyzed:")
    print(input.sample.data.df)
    
    ## Parse formula
    char.formula.elements <- str_split(as.character(dexseq.exp.design)[2], " \\+ ")[[1]]
    formula.reduced.model <- formula(
        paste(
            "~",
            paste(char.formula.elements[1:length(char.formula.elements) - 1],
                  collapse  = " + ")
        ))

    print(paste("Analyzing", sample.base.name))
    print(
        paste(
            "With the full model formula of",
            Reduce(paste, deparse(dexseq.exp.design))
        ))
    print(
        paste(
            "With the reduced model formula of",
            Reduce(paste, deparse(formula.reduced.model))
        ))
    
    ## Apply DRIMSeq filtering as in PMID: 30356428
    count.df <- cbind(
        data.frame(
            gene_id = str_split_fixed(count.dt[, get(ref.column.name)], "_", n = 2)[, 1],
            feature_id = count.dt[, get(ref.column.name)]
        ),
        count.dt[, rownames(input.sample.data.df), with = FALSE]
    )

    replicate.n <- input.sample.data.df[input.sample.data.df$fraction == "ribo1", ] %>%
        nrow
    min.rep <- (replicate.n / 2) %>%
        round %>%
        {min(., 3)}

    ## min.rep <- 6
    
    drimseq.count.data <- dmDSdata(
        counts = count.df,
        samples = data.frame(sample_id = rownames(input.sample.data.df)) %>%
            cbind(., input.sample.data.df)
    ) %>%
        dmFilter(
            min_samps_feature_expr = min.rep,
            min_feature_expr = 10, # TSS has at least 10 count in min.rep samples
            min_samps_feature_prop = min.rep,
            min_feature_prop = 0.05, # TSS expression ratio is at least 0.05% in min.rep samples
            min_samps_gene_expr = min.rep,
            min_gene_expr = 10 # Gene should have count 10 or more in min.rep samples
        )

    dxd <- DEXSeqDataSet(
        countData = as.matrix(counts(drimseq.count.data)[,-c(1:2)]),
        sampleData = input.sample.data.df,
        design = dexseq.exp.design,
        groupID = counts(drimseq.count.data)$gene_id,
        featureID = counts(drimseq.count.data)$feature_id
    ) %>%
        estimateSizeFactors %>%
    estimateDispersions(BPPARAM = BPPARAM, quiet = TRUE) %>%
    testForDEU(
        reducedModel = formula.reduced.model,
        BPPARAM = BPPARAM
    )
        
    dxd.res <- DEXSeqResults(dxd)
    res.dt <- data.table(as.data.frame(dxd.res))

    setnames(
        res.dt,
        old = c("groupID", "featureID"),
        new = c("gene_id", ref.column.name)
    )

    res.dt <- res.dt[, c(
        ref.column.name, "pvalue", "padj"
    ), with = FALSE]

    ## Control gene-wise and tx-wise FDR
    per.gene.qval <- perGeneQValue(dxd.res)
    tx.p.mat <- matrix(dxd.res$pvalue, ncol = 1)
    dimnames(tx.p.mat) <- list(dxd.res$featureID, "transcript")

    tx2gene.df <- as.data.frame(dxd.res[, c("featureID", "groupID")])

    stageRObj <- stageRTx(
        pScreen = per.gene.qval,
        pConfirmation = tx.p.mat,
        pScreenAdjusted = TRUE,
        tx2gene = tx2gene.df
    ) %>%
        stageWiseAdjustment(
            method = "dtu",
            alpha = 0.1
        )

    gene.tx.padj.dt <- data.table(getAdjustedPValues(
        stageRObj,
        order = FALSE,
        onlySignificantGenes = FALSE
    ))

    setnames(
        gene.tx.padj.dt,
        old = c("geneID", "txID", "gene", "transcript"),
        new = c("gene_id", ref.column.name, "gene_FDR", "tx_FDR")
    )

    res.dt <- merge(
        gene.tx.padj.dt,
        res.dt,
        by = ref.column.name
    )
    
    if(nrow(annot.dt) != 0){
        res4export.dt <- merge(
            annot.dt,
            res.dt,
            by = c(ref.column.name, "gene_id")
        )
    } else {
        res4export.dt <- res.dt
    }

    dte.isoform.for.export.file <- file.path(
        s8.2.dte.iso.dir,
        paste0(
            gsub("_NA_\\[\\[\\:digit\\:\\]\\]", "", sample.base.name) %>%
           {gsub("\\(\\[\\[\\:digit\\:\\]\\]_NA\\|\\(3\\|4\\)_Torin1\\)", "__NAvsTorin1", .)},
            ".csv"
        )
    )
    
    fwrite(res4export.dt, file = dte.isoform.for.export.file)
    
    return()
}
```

# Identification of differentially translated mRNA isoforms within a gene

## RCC4 VHL

``` r
sample.base.name <- "RCC4_VHL_EIF4E2_NA_[[:digit:]]_NA"

temp <- analyzeDtIsoforms(
    sample.base.name = sample.base.name,
    count.dt = poly.count.per.tss.dt,
    annot.dt = annot.dt,
    input.sample.data.df = poly.coldata.df %>%
        .[grepl(sample.base.name, rownames(poly.coldata.df)), ],
    dexseq.exp.design = ~ sample + exon + fraction:exon,
    s8.2.dte.iso.dir = s8.2.dte.iso.dir,
    ref.column.name = "tss_name",
    BPPARAM = BPPARAM
)
```

    ## [1] "List of samples analyzed:"
    ##                               cell VHL EIF4E2 gRNA_id clone treatment fraction
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo1 RCC4 VHL EIF4E2      NA    11        NA    ribo1
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo2 RCC4 VHL EIF4E2      NA    11        NA    ribo2
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo3 RCC4 VHL EIF4E2      NA    11        NA    ribo3
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo4 RCC4 VHL EIF4E2      NA    11        NA    ribo4
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo5 RCC4 VHL EIF4E2      NA    11        NA    ribo5
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo6 RCC4 VHL EIF4E2      NA    11        NA    ribo6
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo7 RCC4 VHL EIF4E2      NA    11        NA    ribo7
    ## RCC4_VHL_EIF4E2_NA_1_NA_ribo8 RCC4 VHL EIF4E2      NA    11        NA    ribo8
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo1 RCC4 VHL EIF4E2      NA    13        NA    ribo1
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo2 RCC4 VHL EIF4E2      NA    13        NA    ribo2
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo3 RCC4 VHL EIF4E2      NA    13        NA    ribo3
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo4 RCC4 VHL EIF4E2      NA    13        NA    ribo4
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo5 RCC4 VHL EIF4E2      NA    13        NA    ribo5
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo6 RCC4 VHL EIF4E2      NA    13        NA    ribo6
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo7 RCC4 VHL EIF4E2      NA    13        NA    ribo7
    ## RCC4_VHL_EIF4E2_NA_3_NA_ribo8 RCC4 VHL EIF4E2      NA    13        NA    ribo8
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo1 RCC4 VHL EIF4E2      NA    14        NA    ribo1
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo2 RCC4 VHL EIF4E2      NA    14        NA    ribo2
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo3 RCC4 VHL EIF4E2      NA    14        NA    ribo3
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo4 RCC4 VHL EIF4E2      NA    14        NA    ribo4
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo5 RCC4 VHL EIF4E2      NA    14        NA    ribo5
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo6 RCC4 VHL EIF4E2      NA    14        NA    ribo6
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo7 RCC4 VHL EIF4E2      NA    14        NA    ribo7
    ## RCC4_VHL_EIF4E2_NA_4_NA_ribo8 RCC4 VHL EIF4E2      NA    14        NA    ribo8
    ## [1] "Analyzing RCC4_VHL_EIF4E2_NA_[[:digit:]]_NA"
    ## [1] "With the full model formula of ~sample + exon + fraction:exon"
    ## [1] "With the reduced model formula of ~sample + exon"

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(rse, design, ignoreRank = TRUE): some variables in
    ## design formula are characters, converting to factors

    ## The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of 10%. If a different target OFDR level is of interest,the entire adjustment should be re-run.

## RCC4 noVHL

``` r
sample.base.name <- "RCC4_noVHL_EIF4E2_NA_[[:digit:]]_NA"

temp <- analyzeDtIsoforms(
    sample.base.name = sample.base.name,
    count.dt = poly.count.per.tss.dt,
    annot.dt = annot.dt,
    input.sample.data.df = poly.coldata.df %>%
        .[grepl(sample.base.name, rownames(poly.coldata.df)), ],
    dexseq.exp.design = ~ sample + exon + fraction:exon,
    s8.2.dte.iso.dir = s8.2.dte.iso.dir,
    ref.column.name = "tss_name",
    BPPARAM = BPPARAM
)
```

    ## [1] "List of samples analyzed:"
    ##                                 cell   VHL EIF4E2 gRNA_id clone treatment
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo1 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo2 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo3 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo4 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo5 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo6 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo7 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo8 RCC4 noVHL EIF4E2      NA     1        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo1 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo2 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo3 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo4 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo5 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo6 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo7 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo8 RCC4 noVHL EIF4E2      NA     3        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo1 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo2 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo3 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo4 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo5 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo6 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo7 RCC4 noVHL EIF4E2      NA     4        NA
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo8 RCC4 noVHL EIF4E2      NA     4        NA
    ##                                 fraction
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo1    ribo1
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo2    ribo2
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo3    ribo3
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo4    ribo4
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo5    ribo5
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo6    ribo6
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo7    ribo7
    ## RCC4_noVHL_EIF4E2_NA_1_NA_ribo8    ribo8
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo1    ribo1
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo2    ribo2
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo3    ribo3
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo4    ribo4
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo5    ribo5
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo6    ribo6
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo7    ribo7
    ## RCC4_noVHL_EIF4E2_NA_3_NA_ribo8    ribo8
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo1    ribo1
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo2    ribo2
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo3    ribo3
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo4    ribo4
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo5    ribo5
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo6    ribo6
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo7    ribo7
    ## RCC4_noVHL_EIF4E2_NA_4_NA_ribo8    ribo8
    ## [1] "Analyzing RCC4_noVHL_EIF4E2_NA_[[:digit:]]_NA"
    ## [1] "With the full model formula of ~sample + exon + fraction:exon"
    ## [1] "With the reduced model formula of ~sample + exon"

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(rse, design, ignoreRank = TRUE): some variables in
    ## design formula are characters, converting to factors

    ## The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of 10%. If a different target OFDR level is of interest,the entire adjustment should be re-run.

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
    ##  [1] knitr_1.28                  stringr_1.4.0              
    ##  [3] magrittr_1.5                data.table_1.12.8          
    ##  [5] dplyr_1.0.0                 khroma_1.3.0               
    ##  [7] ggplot2_3.3.1               DEXSeq_1.34.0              
    ##  [9] RColorBrewer_1.1-2          AnnotationDbi_1.50.0       
    ## [11] DESeq2_1.28.0               BiocParallel_1.22.0        
    ## [13] stageR_1.10.0               SummarizedExperiment_1.18.1
    ## [15] DelayedArray_0.14.0         matrixStats_0.56.0         
    ## [17] Biobase_2.48.0              GenomicRanges_1.40.0       
    ## [19] GenomeInfoDb_1.24.0         IRanges_2.22.1             
    ## [21] S4Vectors_0.26.0            BiocGenerics_0.34.0        
    ## [23] DRIMSeq_1.16.0              rmarkdown_2.2              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2             edgeR_3.30.0           bit64_0.9-7           
    ##  [4] splines_4.0.0          assertthat_0.2.1       askpass_1.1           
    ##  [7] statmod_1.4.34         BiocFileCache_1.12.0   blob_1.2.1            
    ## [10] Rsamtools_2.4.0        GenomeInfoDbData_1.2.3 yaml_2.2.1            
    ## [13] progress_1.2.2         pillar_1.4.4           RSQLite_2.2.0         
    ## [16] lattice_0.20-41        glue_1.4.1             limma_3.44.1          
    ## [19] digest_0.6.25          XVector_0.28.0         colorspace_1.4-1      
    ## [22] htmltools_0.4.0        Matrix_1.2-18          plyr_1.8.6            
    ## [25] XML_3.99-0.3           pkgconfig_2.0.3        biomaRt_2.44.0        
    ## [28] genefilter_1.70.0      zlibbioc_1.34.0        purrr_0.3.4           
    ## [31] xtable_1.8-4           scales_1.1.1           openssl_1.4.1         
    ## [34] tibble_3.0.1           annotate_1.66.0        generics_0.0.2        
    ## [37] ellipsis_0.3.1         withr_2.4.1            survival_3.1-12       
    ## [40] crayon_1.3.4           memoise_1.1.0          evaluate_0.14         
    ## [43] hwriter_1.3.2          tools_4.0.0            prettyunits_1.1.1     
    ## [46] hms_0.5.3              lifecycle_0.2.0        munsell_0.5.0         
    ## [49] locfit_1.5-9.4         Biostrings_2.56.0      compiler_4.0.0        
    ## [52] rlang_0.4.10           grid_4.0.0             RCurl_1.98-1.2        
    ## [55] rappdirs_0.3.1         bitops_1.0-6           gtable_0.3.0          
    ## [58] curl_4.3               DBI_1.1.0              reshape2_1.4.4        
    ## [61] R6_2.4.1               bit_1.1-15.2           stringi_1.4.6         
    ## [64] Rcpp_1.0.4.6           vctrs_0.3.1            geneplotter_1.66.0    
    ## [67] dbplyr_1.4.4           tidyselect_1.1.0       xfun_0.14
