---
title: "s7-2 Define HIF1A and HIF2A binding ratio"
author: "Yoichiro Sugimoto"
date: "13 November, 2021"
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

Based on ChIP-Seq data, HIF binding sites and the ratio of HIF2A binding relative to that of HIF1A will be defined.


# Environment setup and data preprocessing


```r
## Bioconductor
library("GenomicFeatures")
library("DiffBind")
## rmarkdown
library("knitr")
library("kableExtra")

## Specify the number of CPUs to be used
processors <- 8

## library("BiocParallel")
## register(MulticoreParam(processors))

temp <- sapply(list.files("../functions", full.names = TRUE), source)

results.dir <- file.path("../../results")

s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")
s4.3.tx.info.dir <- file.path(s4.tss.dir, "s4-3-transcript-info")
s4.3.1.tx.info.rcc4.dir <- file.path(s4.3.tx.info.dir, "s4-3-1-transcript-info-for-RCC4")

s7.dir <- file.path(results.dir, "s7-HIF-binding-site")
s7.1.dir <- file.path(s7.dir, "s7-1-ChIP-Seq")
s7.1.macs.dir <- file.path(s7.1.dir, "macs_out")
s7.1.hre.bed.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_HRE_bed")
s7.1.hre.narrowPeaks.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_HRE_narrowPeaks")
s7.2.hif1.2a.dir <- file.path(s7.dir, "s7-2-HIF1A-and-HIF2A-ratio")

create.dirs(c(
    s7.2.hif1.2a.dir
))

set.seed(0)
```


# Define HIF binding sites and the ratio of HIF2A/HIF1A binding



```r
reprocessed.chip.seq.data.dir <- file.path(
    "../../../20200519_reprocessing_HIF_ChIP-Seq"
)

chip.seq.sample.dt <- file.path(
    "../../data/sample_data/20210203_ChIPSeq_sample_data.csv"
) %>% fread

chip.seq.sample.dt[, `:=`(
    bamReads = file.path(
        reprocessed.chip.seq.data.dir, "processed_data/sorted_bam", bamReads
    ),
    Peaks = file.path(
        s7.1.hre.narrowPeaks.dir, Peaks
    )
)]

all.chip.peak.gr <- dba(sampleSheet = chip.seq.sample.dt) %>%
    dba.count %>%
    dba.contrast(categories = DBA_FACTOR, minMembers = 2) %>%
    dba.analyze(bFullLibrarySize = FALSE) %>%
    dba.report(th = 1)
```

```
## HIF1A_1  HIF1A   1 narrow
```

```
## HIF1A_2  HIF1A   2 narrow
```

```
## HIF2A_1  HIF2A   1 narrow
```

```
## HIF2A_2  HIF2A   2 narrow
```

```
## converting counts to integer mode
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```r
readChIPSeqData <- function(chip.seq.core.name){
    idr.narrowpeak.file <- file.path(
        s7.1.hre.narrowPeaks.dir,
        paste0(chip.seq.core.name, ".narrowPeak")
    )
    temp.idr.narrowpeak.dt <- fread(idr.narrowpeak.file)
    
    idr.narrowpeak.dt <- data.table(
        chr = temp.idr.narrowpeak.dt[, V1],
        start = temp.idr.narrowpeak.dt[, V2 + V10] + 1, # bed is 0-base
        end = temp.idr.narrowpeak.dt[, V2 + V10], 
        enrichment = temp.idr.narrowpeak.dt[, V7],
        mlog10qval = temp.idr.narrowpeak.dt[, V9],
        peak = temp.idr.narrowpeak.dt[, V2 + V10] + 1 # bed is 0-base
    )

    chip.seq.summit.gr <- makeGRangesFromDataFrame(
        idr.narrowpeak.dt, keep.extra.columns = TRUE
    )
}

hif1a.chip.seq.summit.gr <- readChIPSeqData("RCC4_N_HIF1A")
hif2a.chip.seq.summit.gr <- readChIPSeqData("RCC4_N_HIF2A")

hif1a.ol.dt <- findOverlaps(
    hif1a.chip.seq.summit.gr,
    all.chip.peak.gr,
) %>% as.data.frame %>% data.table

hif2a.ol.dt <- findOverlaps(
    hif2a.chip.seq.summit.gr,
    all.chip.peak.gr
) %>% as.data.frame %>% data.table

mcols(all.chip.peak.gr)[["index"]] <- 1:length(all.chip.peak.gr)

all.chip.peak.mcols.dt <- all.chip.peak.gr %>%
    as.data.frame %>% data.table

np.cols <- c("enrichment", "mlog10qval", "peak")

h1.mcol.dt <- cbind(
    all.chip.peak.mcols.dt[hif1a.ol.dt[, subjectHits], "index", with = FALSE],
    mcols(hif1a.chip.seq.summit.gr[hif1a.ol.dt[, queryHits]]) %>% as.data.frame %>% data.table
)
setnames(h1.mcol.dt, old = np.cols, new = paste0("HIF1A_", np.cols))

h2.mcol.dt <- cbind(
    all.chip.peak.mcols.dt[hif2a.ol.dt[, subjectHits], "index", with = FALSE],
    mcols(hif2a.chip.seq.summit.gr[hif2a.ol.dt[, queryHits]]) %>% as.data.frame %>% data.table
)
setnames(h2.mcol.dt, old = np.cols, new = paste0("HIF2A_", np.cols))

all.chip.peak.comp.mcols.dt <- Reduce(
    function(...) merge(..., all = TRUE, by = "index"),
    list(
        all.chip.peak.mcols.dt,
        h1.mcol.dt,
        h2.mcol.dt
    )) %>%
    {.[order(index, HIF1A_enrichment, HIF2A_enrichment)][!duplicated(index)]}

all.chip.peak.comp.mcols.dt[, `:=`(
    peak_position = case_when(
        is.na(HIF2A_peak) & !is.na(HIF1A_peak) ~ HIF1A_peak,
        is.na(HIF1A_peak) & !is.na(HIF2A_peak) ~ HIF2A_peak,
        ## HIF1A_mlog10qval >= HIF2A_mlog10qval ~ HIF1A_peak,
        Fold >= 0 ~ HIF1A_peak,
        Fold < 0 ~ HIF2A_peak
    )
)]

fwrite(
    all.chip.peak.comp.mcols.dt,
    file.path(
        s7.2.hif1.2a.dir,
        "HIF1A-HIF2A-binding-summary.csv"
    )
)
```

## Export HIF binding sites in bed format



```r
hif.bed.dt <- all.chip.peak.comp.mcols.dt[, list(
    chrom = seqnames,
    start = start - 1,
    end = end,
    name = paste0("peak_", index, "_HIF2A/HIF1_log2fc_", - Fold),
    score = 1,
    strand = strand,
    thick_start = peak_position - 1,
    thick_end = peak_position,
    itemRgb = "0,0,255",
    block_count = 1,
    blockSizes = end - start - 1,
    blockStart = 0    
)]

fwrite(
    hif.bed.dt,
    file.path(
        s7.2.hif1.2a.dir,
        "HIF-binding-site.bed"
    ),
    col.names = FALSE, sep = "\t"
)
```



# Assign closest HIF binding sites to all TSS



```r
filtered.tss.with.quantile.dt <- file.path(
  s4.1.6.filtered.tss.dir,
  "filtered-tss-with-quantile.csv"
) %>% fread

filtered.tss.with.quantile.dt <- filtered.tss.with.quantile.dt[gene_id != ""]

filtered.tss.with.quantile.gr <- makeGRangesFromDataFrame(
    filtered.tss.with.quantile.dt,
    keep.extra.columns = TRUE
)

chip.seq.summit.gr <- makeGRangesFromDataFrame(
    all.chip.peak.comp.mcols.dt[!is.na(peak_position)],
    keep.extra.columns = TRUE
)

start(chip.seq.summit.gr) <- mcols(chip.seq.summit.gr)[["peak_position"]]
end(chip.seq.summit.gr) <- mcols(chip.seq.summit.gr)[["peak_position"]] 

nearest.hif.dt <- nearest(
    filtered.tss.with.quantile.gr,
    chip.seq.summit.gr,
    ignore.strand = TRUE,
    select = "all"
) %>% as.data.frame %>% data.table

tss.hif.pos.dt <- cbind(
    filtered.tss.with.quantile.dt[nearest.hif.dt[, queryHits]],
    as.data.frame(chip.seq.summit.gr[nearest.hif.dt[, subjectHits]]) %>%
    data.table %>% {
        .[, list(
            hre_position = start,
            Conc_HIF1A, Conc_HIF2A, Fold
        )]
    }
)

tss.hif.pos.dt[, `:=`(
    tss_pos_to_hif = ifelse(strand == "+", 1 , -1) * (q50 - hre_position)
)]
tss.hif.pos.dt[, `:=`(
    log10_tss_pos_to_hif = sign(tss_pos_to_hif) * log10(abs(tss_pos_to_hif))
)]
tss.hif.pos.dt <- tss.hif.pos.dt[gene_id != ""]

tss.hif.pos.dt[, `:=`(
    HIF2A_enrichment_with_dist_th = if_else(
        abs(tss_pos_to_hif) < 500 * 10^3, - Fold, NA_real_
    )
)]

fwrite(
    tss.hif.pos.dt,
    file.path(
        s7.2.hif1.2a.dir,
        "filtered-tss-and-nearest-hif-binding-position.csv"
    )
)
```


# session information


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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] stringr_1.4.0               magrittr_1.5               
##  [3] data.table_1.12.8           dplyr_1.0.0                
##  [5] khroma_1.3.0                ggplot2_3.3.1              
##  [7] kableExtra_1.1.0            knitr_1.28                 
##  [9] DiffBind_2.16.0             SummarizedExperiment_1.18.1
## [11] DelayedArray_0.14.0         matrixStats_0.56.0         
## [13] GenomicFeatures_1.40.0      AnnotationDbi_1.50.0       
## [15] Biobase_2.48.0              GenomicRanges_1.40.0       
## [17] GenomeInfoDb_1.24.0         IRanges_2.22.1             
## [19] S4Vectors_0.26.0            BiocGenerics_0.34.0        
## [21] rmarkdown_2.2              
## 
## loaded via a namespace (and not attached):
##   [1] amap_0.8-18              colorspace_1.4-1         rjson_0.2.20            
##   [4] hwriter_1.3.2            ellipsis_0.3.1           XVector_0.28.0          
##   [7] rstudioapi_0.11          ggrepel_0.8.2            bit64_0.9-7             
##  [10] xml2_1.3.2               splines_4.0.0            geneplotter_1.66.0      
##  [13] jsonlite_1.7.2           Rsamtools_2.4.0          annotate_1.66.0         
##  [16] GO.db_3.11.4             dbplyr_1.4.4             png_0.1-7               
##  [19] pheatmap_1.0.12          graph_1.66.0             readr_1.3.1             
##  [22] compiler_4.0.0           httr_1.4.2               GOstats_2.54.0          
##  [25] backports_1.1.7          assertthat_0.2.1         Matrix_1.2-18           
##  [28] limma_3.44.1             htmltools_0.4.0          prettyunits_1.1.1       
##  [31] tools_4.0.0              gtable_0.3.0             glue_1.4.1              
##  [34] GenomeInfoDbData_1.2.3   Category_2.54.0          systemPipeR_1.22.0      
##  [37] rsvg_2.1                 batchtools_0.9.14        rappdirs_0.3.1          
##  [40] V8_3.2.0                 ShortRead_1.46.0         Rcpp_1.0.4.6            
##  [43] vctrs_0.3.1              Biostrings_2.56.0        rtracklayer_1.48.0      
##  [46] xfun_0.14                rvest_0.3.5              lifecycle_0.2.0         
##  [49] gtools_3.8.2             XML_3.99-0.3             edgeR_3.30.0            
##  [52] zlibbioc_1.34.0          scales_1.1.1             BSgenome_1.56.0         
##  [55] VariantAnnotation_1.34.0 hms_0.5.3                RBGL_1.64.0             
##  [58] RColorBrewer_1.1-2       yaml_2.2.1               curl_4.3                
##  [61] memoise_1.1.0            biomaRt_2.44.0           latticeExtra_0.6-29     
##  [64] stringi_1.4.6            RSQLite_2.2.0            genefilter_1.70.0       
##  [67] checkmate_2.0.0          caTools_1.18.0           BiocParallel_1.22.0     
##  [70] DOT_0.1                  rlang_0.4.10             pkgconfig_2.0.3         
##  [73] bitops_1.0-6             evaluate_0.14            lattice_0.20-41         
##  [76] purrr_0.3.4              GenomicAlignments_1.24.0 bit_1.1-15.2            
##  [79] tidyselect_1.1.0         GSEABase_1.50.0          AnnotationForge_1.30.1  
##  [82] DESeq2_1.28.0            R6_2.4.1                 gplots_3.1.1            
##  [85] generics_0.0.2           base64url_1.4            DBI_1.1.0               
##  [88] pillar_1.4.4             withr_2.4.1              survival_3.1-12         
##  [91] RCurl_1.98-1.2           tibble_3.0.1             crayon_1.3.4            
##  [94] KernSmooth_2.23-17       BiocFileCache_1.12.0     jpeg_0.1-8.1            
##  [97] progress_1.2.2           locfit_1.5-9.4           grid_4.0.0              
## [100] blob_1.2.1               Rgraphviz_2.32.0         webshot_0.5.2           
## [103] digest_0.6.25            xtable_1.8-4             brew_1.0-6              
## [106] openssl_1.4.1            munsell_0.5.0            viridisLite_0.3.0       
## [109] askpass_1.1
```
