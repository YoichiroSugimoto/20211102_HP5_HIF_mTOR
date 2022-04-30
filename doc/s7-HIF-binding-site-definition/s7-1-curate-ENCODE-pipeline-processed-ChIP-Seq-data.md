s7-1 Curate ENCODE pipeline processed ChIP-Seq data
================
Yoichiro Sugimoto
29 April, 2022

  - [Overview](#overview)
  - [Environment setup and data
    preprocessing](#environment-setup-and-data-preprocessing)
  - [Covert IDR conservative narrowPeak file to bed
    file](#covert-idr-conservative-narrowpeak-file-to-bed-file)
      - [Generate IDR non duplicated
        peaks](#generate-idr-non-duplicated-peaks)
  - [Rename fold change bigwig file](#rename-fold-change-bigwig-file)
  - [Peak calling with MACS2](#peak-calling-with-macs2)
      - [Run MACS2 for pooled data](#run-macs2-for-pooled-data)
      - [Filter MACS output with IDR](#filter-macs-output-with-idr)
  - [Identify closest HRE to the cluster
    centre](#identify-closest-hre-to-the-cluster-centre)
  - [Session information](#session-information)

# Overview

The intersections of ENCODE ChIP-Seq pipeline identified peaks and MACS
identified peaks will be used for the following analyses.

# Environment setup and data preprocessing

``` r
library("GenomicRanges")

## Specify the number of CPUs to be used
processors <- 8

temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)

set.seed(0)
```

``` r
results.dir <- file.path("../../results")

s7.dir <- file.path(results.dir, "s7-HIF-binding-site")
s7.1.dir <- file.path(s7.dir, "s7-1-ChIP-Seq")

s7.1.peak.bed.dir <- file.path(s7.1.dir, "peak_bed")
s7.1.idr.filtered.narrowPeak.dir <- file.path(s7.1.dir, "IDR_filtered_narrowPeaks")
s7.1.idr.filtered.bed.dir <- file.path(s7.1.dir, "IDR_filtered_bed")
s7.1.fc.bigwig.dir <- file.path(s7.1.dir, "foldchange_bigwig")
s7.1.macs.dir <- file.path(s7.1.dir, "macs_out")
s7.1.macs.bed.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_bed")
s7.1.macs.narrowPeaks.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_narrowPeaks")
s7.1.hre.bed.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_HRE_bed")
s7.1.hre.narrowPeaks.dir <- file.path(s7.1.macs.dir, "IDR_filtered_macs_HRE_narrowPeaks")


create.dirs(c(
    s7.dir,
    s7.1.dir,
    s7.1.peak.bed.dir,
    s7.1.fc.bigwig.dir,
    s7.1.macs.dir,
    s7.1.idr.filtered.narrowPeak.dir,
    s7.1.idr.filtered.bed.dir,
    s7.1.macs.bed.dir,
    s7.1.macs.narrowPeaks.dir,
    s7.1.hre.bed.dir,
    s7.1.hre.narrowPeaks.dir 
))
```

# Covert IDR conservative narrowPeak file to bed file

``` r
chip.encodepl.dir <- file.path("../../../20200519_reprocessing_HIF_ChIP-Seq")

covertNp2Bed <- function(chip.seq.exp, chip.encodepl.dir){

    idr.narrowpeak.file <- file.path(
        list.dirs(file.path(chip.encodepl.dir, chip.seq.exp, "chip"), recursive = FALSE),
        "call-reproducibility_idr/execution/idr.conservative_peak.regionPeak.gz"
    )
    temp.idr.narrowpeak.dt <- fread(cmd = paste("zcat", idr.narrowpeak.file))

    idr.bed.dt <- cbind(
        temp.idr.narrowpeak.dt[, 1:6],
        data.table(
            thick_start = temp.idr.narrowpeak.dt[, V2 + V10],
            thick_end = temp.idr.narrowpeak.dt[, V2 + V10 + 1],
            itemRgb = "0,0,255",
            block_count = 1,
            blockSizes = temp.idr.narrowpeak.dt[, V3 - V2],
            blockStart = 0
        )
    )

    idr.bed.dt <- idr.bed.dt[order(V2)][
        order(match(V1, paste0("chr", c(1:22, "X", "Y"))))
    ]
    
    out.peak.file.name <- file.path(
        s7.1.peak.bed.dir,
        paste0(chip.seq.exp, ".bed")
    )

    fwrite(idr.bed.dt, file = out.peak.file.name, col.names = FALSE, sep = "\t")

    return()

}

temp <- lapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"),
    covertNp2Bed,
    chip.encodepl.dir = chip.encodepl.dir
)
```

## Generate IDR non duplicated peaks

Some of IDR identified peaks are overlapping each other. Here I filter
out such peaks.

``` r
filterIDRselectedPeaks <- function(chip.seq.exp, chip.encodepl.dir){

    idr.narrowpeak.file <- file.path(
        list.dirs(file.path(chip.encodepl.dir, chip.seq.exp, "chip"), recursive = FALSE),
        "call-reproducibility_idr/execution/idr.conservative_peak.regionPeak.gz"
    )
    temp.idr.narrowpeak.dt <- fread(cmd = paste("zcat", idr.narrowpeak.file))

    setnames(
        temp.idr.narrowpeak.dt,
        old = paste0("V", 1:10),
        new = c("chr", "peak_start", "peak_end", "peak_name", "score", "strand",
                "fold_change", "mlog10p", "mlog10q", "rel_summit_pos")
    )
    temp.idr.narrowpeak.dt[, summit_pos := peak_start + rel_summit_pos]
    
    temp.idr.narrowpeak.gr <- makeGRangesFromDataFrame(
        temp.idr.narrowpeak.dt,
        start.field = "peak_start",
        end.field = "peak_end",
        keep.extra.columns = TRUE
    )

    ol.dt <- data.table(as.data.frame(findOverlaps(
        temp.idr.narrowpeak.gr,
        reduce(temp.idr.narrowpeak.gr)
    )))

    temp.idr.narrowpeak.dt <- cbind(
        temp.idr.narrowpeak.dt[ol.dt[, queryHits]],
        data.table(rd_idx = ol.dt[, subjectHits])
    )

    temp.idr.narrowpeak.dt <- temp.idr.narrowpeak.dt[
        temp.idr.narrowpeak.dt[, .I[fold_change == max(fold_change)], by = rd_idx]$V1
    ]
    
    if(nrow(temp.idr.narrowpeak.dt[duplicated(rd_idx)]) != 0){
        stop("Multiple summits within a reduced IDR peak")
    } else {}

    temp.idr.narrowpeak.dt[, `:=`(
        summit_pos = NULL,
        rd_idx = NULL
    )]

    narrowPeak.out <- file.path(
        s7.1.idr.filtered.narrowPeak.dir,
        paste0(chip.seq.exp, ".narrowPeak")
    )

    fwrite(temp.idr.narrowpeak.dt, file = narrowPeak.out, sep = "\t", col.names = FALSE)

    ## Convert to bed
    idr.bed.dt <- cbind(
        temp.idr.narrowpeak.dt[, 1:6],
        data.table(
            thick_start = temp.idr.narrowpeak.dt[, peak_start + rel_summit_pos],
            thick_end = temp.idr.narrowpeak.dt[, peak_start + rel_summit_pos + 1],
            itemRgb = "0,0,255",
            block_count = 1,
            blockSizes = temp.idr.narrowpeak.dt[, peak_end - peak_start],
            blockStart = 0
        )
    )

    idr.bed.dt <- idr.bed.dt[order(peak_start)][
        order(match(chr, paste0("chr", c(1:22, "X", "Y"))))
    ]
    
    out.peak.file.name <- file.path(
        s7.1.idr.filtered.bed.dir,
        paste0(chip.seq.exp, ".bed")
    )

    fwrite(idr.bed.dt, file = out.peak.file.name, col.names = FALSE, sep = "\t")

    return()

}

temp <- lapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"),
    filterIDRselectedPeaks,
    chip.encodepl.dir = chip.encodepl.dir
)
```

# Rename fold change bigwig file

``` r
renameBw <- function(chip.seq.exp, chip.encodepl.dir){
    pooled.fc.bw.file <- file.path(
        list.dirs(file.path(chip.encodepl.dir, chip.seq.exp, "chip"), recursive = FALSE),
        "call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.fc.signal.bigwig"
    )

    out.bw.file <- file.path(
        s7.1.fc.bigwig.dir,
        paste0(chip.seq.exp, ".fc.bigwig")
    )

    file.copy(from = pooled.fc.bw.file, to = out.bw.file, overwrite = TRUE)
}

temp <- mclapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"),
    renameBw,
    chip.encodepl.dir = chip.encodepl.dir,
    mc.cores = processors
)
```

# Peak calling with MACS2

IDR peak sometimes contain multiple peaks in a peak. Thus peak will be
called for pooled sample to be compared.

## Run MACS2 for pooled data

``` r
runMACS2 <- function(chip.seq.exp, chip.encodepl.dir, s7.1.macs.dir){

    temp.bam.files <- list.files(
        file.path(
            list.dirs(file.path(chip.encodepl.dir, chip.seq.exp, "chip"), recursive = FALSE),
            "call-filter"
        ),
        pattern = "nodup.bam$", recursive = TRUE, full.names = TRUE
    )
    bam.files <- temp.bam.files[!grepl("glob", temp.bam.files)]

    temp.ctrl.bam.files <- list.files(
        file.path(
            list.dirs(file.path(chip.encodepl.dir, chip.seq.exp, "chip"), recursive = FALSE),
            "call-filter_ctl"
        ),
        pattern = "nodup.bam$", recursive = TRUE, full.names = TRUE
    )
    ctrl.bam.files <- temp.ctrl.bam.files[!grepl("glob", temp.ctrl.bam.files)]

    macs.cmd <- paste(
        "macs2", "callpeak",
        "-t", paste(bam.files, collapse = " "),
        "-c", paste(ctrl.bam.files, collapse = " "),
        "-f", "BAMPE",
        "-g", "hs",
        "-n", chip.seq.exp,
        "-q", 0.1,
        "--call-summits",
        "--outdir", file.path(s7.1.macs.dir, chip.seq.exp)
    )

    system.cat(macs.cmd)
    
    return()

}

temp <- mclapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"),
    runMACS2,
    chip.encodepl.dir = chip.encodepl.dir,
    s7.1.macs.dir = s7.1.macs.dir,
    mc.cores = processors
)
```

## Filter MACS output with IDR

``` r
filterMACS2out <- function(chip.seq.exp, chip.encodepl.dir, s7.1.macs.dir){

    macs2.out <- file.path(
        s7.1.macs.dir, chip.seq.exp,
        paste0(chip.seq.exp, "_peaks.narrowPeak")
        )

    idr.out <- file.path(s7.1.peak.bed.dir, paste0(chip.seq.exp, ".bed"))

    macs2.dt <- fread(macs2.out)
    idr.dt <- fread(idr.out)

    setnames(
        macs2.dt,
        old = paste0("V", 1:10),
        new = c("chr", "peak_start", "peak_end", "peak_name", "score", "strand",
                "fold_change", "mlog10p", "mlog10q", "rel_summit_pos")
    )
    macs2.dt[, summit_pos := peak_start + rel_summit_pos]

    macs2.gr <- makeGRangesFromDataFrame(
        macs2.dt, start.field = "summit_pos", end.field = "summit_pos",
        keep.extra.columns = TRUE
    )

    setnames(
        idr.dt,
        old = paste0("V", 1:3),
        new = c("chr", "start", "end")
    )

    idr.gr <- makeGRangesFromDataFrame(idr.dt)
    rd.idr.gr <- reduce(idr.gr)
    
    ol.dt <- data.table(
        as.data.frame(findOverlaps(
            macs2.gr, rd.idr.gr
        ))
    )
    
    idr.filtered.macs2.dt <- cbind(
        macs2.dt[ol.dt[, queryHits]],
        data.table(idr_peak_idx = ol.dt[, subjectHits])
    )

    idr.filtered.macs2.dt <- idr.filtered.macs2.dt[
        idr.filtered.macs2.dt[, .I[mlog10q == max(mlog10q)], by = idr_peak_idx]$V1
    ]

    if(nrow(idr.filtered.macs2.dt[duplicated(idr_peak_idx)]) != 0){
        stop("Multiple summits within a reduced IDR peak")
    } else {}

    idr.filtered.macs2.dt[, `:=`(
        summit_pos = NULL,
        idr_peak_idx = NULL
    )]

    narrowPeak.out <- file.path(
        s7.1.macs.narrowPeaks.dir,
        paste0(chip.seq.exp, ".narrowPeak")
    )

    fwrite(idr.filtered.macs2.dt, file = narrowPeak.out, sep = "\t", col.names = FALSE)

    ## Convert to bed
    idr.bed.dt <- cbind(
        idr.filtered.macs2.dt[, 1:6],
        data.table(
            thick_start = idr.filtered.macs2.dt[, peak_start + rel_summit_pos],
            thick_end = idr.filtered.macs2.dt[, peak_start + rel_summit_pos + 1],
            itemRgb = "0,0,255",
            block_count = 1,
            blockSizes = idr.filtered.macs2.dt[, peak_end - peak_start],
            blockStart = 0
        )
    )

    idr.bed.dt <- idr.bed.dt[order(peak_start)][
        order(match(chr, paste0("chr", c(1:22, "X", "Y"))))
    ]
    
    out.peak.file.name <- file.path(
        s7.1.macs.bed.dir,
        paste0(chip.seq.exp, ".bed")
    )

    fwrite(idr.bed.dt, file = out.peak.file.name, col.names = FALSE, sep = "\t")

    return()

}

temp <- mclapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"), 
    filterMACS2out,
    chip.encodepl.dir = chip.encodepl.dir,
    s7.1.macs.dir = s7.1.macs.dir,
    mc.cores = processors
)
```

# Identify closest HRE to the cluster centre

``` r
library("BSgenome.Hsapiens.UCSC.hg38") 
```

    ## Loading required package: BSgenome

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: rtracklayer

``` r
generatePeakHreBed <- function(chip.seq.exp, s7.1.macs.narrowPeaks.dir, s7.1.hre.narrowPeaks.dir, s7.1.hre.bed.dir){

    chip.peak.dt <- fread(
        file.path(s7.1.macs.narrowPeaks.dir, paste0(chip.seq.exp, ".narrowPeak"))
    )

    setnames(
        chip.peak.dt,
        old = colnames(chip.peak.dt),
        new = c("chr", "start", "end", "peak_id", "score", "strand", "log2fc", "mlog10pvalue", "mlog10qvalue", "rel_summit_position")
    )

    ## Remove overlapping peak
    chip.peak.dt[, nondup_peak_id := stringr::str_extract(peak_id, pattern = paste0(chip.seq.exp, "_peak_[:digit:]+"))]

    chip.peak.dt <- chip.peak.dt[order(mlog10qvalue, log2fc, decreasing = TRUE)][
      , head(.SD, 1), by = nondup_peak_id
    ]

    ## narropeak is bed file and thus 0 base
    hre.search.range <- 50

    chip.peak.dt[, `:=`(
        summit_start = start + rel_summit_position + 1, # Here I use 1 base
        summit_end = start + rel_summit_position + 1
    )] %>%
        {.[, `:=`(
             region_start = summit_start - hre.search.range,
             region_end = summit_end + hre.search.range
         )]}

    chip.summit.gr <- makeGRangesFromDataFrame(
        chip.peak.dt[, .(chr, region_start, region_end, peak_id)],
        start.field = "region_start",
        end.field = "region_end",
        keep.extra.columns = TRUE
    )

    summit.range.seq <- getSeq(
        BSgenome.Hsapiens.UCSC.hg38,
        chip.summit.gr
    )

    names(summit.range.seq) <- mcols(chip.summit.gr)[, "peak_id"]

    hre.motif <- "RCGTG"

    hre.sense.match.dt <- vmatchPattern(
        pattern = hre.motif,
        summit.range.seq,
        fixed = FALSE
    ) %>%
        as.data.frame %>% data.table %>%
    {.[, peak_id := names(summit.range.seq)[group]]}

    hre.rc.match.dt <- vmatchPattern(
        pattern = reverseComplement(DNAString(hre.motif)),
        summit.range.seq,
        fixed = FALSE
    ) %>%
        as.data.frame %>% data.table %>%
    {.[, peak_id := names(summit.range.seq)[group]]}

    hre.match.dt <- rbind(
        hre.sense.match.dt,
        hre.rc.match.dt
    )

    hre.match.dt[, `:=`(
        distance_to_summit = start + 2 - (hre.search.range + 1)
    )]

    closest.hre.positon.dt <- hre.match.dt[order(start)][order(abs(distance_to_summit))][, head(.SD, 1), by = peak_id]

    chip.hre.dt <- merge(
        chip.peak.dt,
        closest.hre.positon.dt[, .(peak_id, distance_to_summit)],
        by = "peak_id"
    )

    chip.hre.dt[, `:=`(
        hre_start = start + rel_summit_position + distance_to_summit, # 0-base
        hre_end = start + rel_summit_position + distance_to_summit + 1
    )]

    chip.hre.narrowpeak.dt <- chip.hre.dt[, .(
        chr, start, end, peak_id, score, strand, log2fc, mlog10pvalue, mlog10qvalue, hre_start - start
    )]

    narrowPeak.out <- file.path(
        s7.1.hre.narrowPeaks.dir,
        paste0(chip.seq.exp, ".narrowPeak")
    )

    fwrite(chip.hre.narrowpeak.dt, file = narrowPeak.out, sep = "\t", col.names = FALSE)

    ## Convert to bed
    hre.bed.dt <- cbind(
        chip.hre.narrowpeak.dt[, 1:6],
        data.table(
            thick_start = chip.hre.narrowpeak.dt[, start + V10],
            thick_end = chip.hre.narrowpeak.dt[, start + V10 + 1],
            itemRgb = "0,0,255",
            block_count = 1,
            blockSizes = chip.hre.narrowpeak.dt[, end - start],
            blockStart = 0
        )
    )

    hre.bed.dt <- hre.bed.dt[order(start)][
        order(match(chr, paste0("chr", c(1:22, "X", "Y"))))
    ]

    out.peak.file.name <- file.path(
        s7.1.hre.bed.dir,
        paste0(chip.seq.exp, ".bed")
    )

    fwrite(hre.bed.dt, file = out.peak.file.name, col.names = FALSE, sep = "\t")

    return()
}


temp <- mclapply(
    c("RCC4_N_HIF1A", "RCC4_N_HIF2A"),
    generatePeakHreBed,
    s7.1.macs.narrowPeaks.dir = s7.1.macs.narrowPeaks.dir,
    s7.1.hre.narrowPeaks.dir = s7.1.hre.narrowPeaks.dir,
    s7.1.hre.bed.dir = s7.1.hre.bed.dir
)
```

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
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.3 BSgenome_1.56.0                  
    ##  [3] rtracklayer_1.48.0                Biostrings_2.56.0                
    ##  [5] XVector_0.28.0                    knitr_1.28                       
    ##  [7] stringr_1.4.0                     magrittr_1.5                     
    ##  [9] data.table_1.12.8                 dplyr_1.0.0                      
    ## [11] khroma_1.3.0                      ggplot2_3.3.1                    
    ## [13] GenomicRanges_1.40.0              GenomeInfoDb_1.24.0              
    ## [15] IRanges_2.22.1                    S4Vectors_0.26.0                 
    ## [17] BiocGenerics_0.34.0               rmarkdown_2.2                    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6                compiler_4.0.0             
    ##  [3] pillar_1.4.4                bitops_1.0-6               
    ##  [5] tools_4.0.0                 zlibbioc_1.34.0            
    ##  [7] digest_0.6.25               lattice_0.20-41            
    ##  [9] evaluate_0.14               lifecycle_0.2.0            
    ## [11] tibble_3.0.1                gtable_0.3.0               
    ## [13] pkgconfig_2.0.3             rlang_0.4.10               
    ## [15] Matrix_1.2-18               DelayedArray_0.14.0        
    ## [17] yaml_2.2.1                  xfun_0.14                  
    ## [19] GenomeInfoDbData_1.2.3      withr_2.4.1                
    ## [21] generics_0.0.2              vctrs_0.3.1                
    ## [23] tidyselect_1.1.0            grid_4.0.0                 
    ## [25] Biobase_2.48.0              glue_1.4.1                 
    ## [27] R6_2.4.1                    BiocParallel_1.22.0        
    ## [29] XML_3.99-0.3                purrr_0.3.4                
    ## [31] matrixStats_0.56.0          GenomicAlignments_1.24.0   
    ## [33] Rsamtools_2.4.0             scales_1.1.1               
    ## [35] htmltools_0.4.0             ellipsis_0.3.1             
    ## [37] SummarizedExperiment_1.18.1 colorspace_1.4-1           
    ## [39] stringi_1.4.6               RCurl_1.98-1.2             
    ## [41] munsell_0.5.0               crayon_1.3.4
