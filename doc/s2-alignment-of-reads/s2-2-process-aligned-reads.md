s2-2 Process aligned reads (UMI accounting and tss count export)
================
Yoichiro Sugimoto
01 March, 2022

  - [Overview](#overview)
  - [Set up](#set-up)
  - [Generate TSS bam file](#generate-tss-bam-file)
  - [Deduplication with UMI](#deduplication-with-umi)
  - [Generate deduplicated bam file](#generate-deduplicated-bam-file)
  - [Create read count table by
    genes](#create-read-count-table-by-genes)
  - [Session Information](#session-information)

# Overview

Process aligned reads to:

  - Generate TSS bam files
  - Account UMIs based on the TSS bam files
  - Generate deduplicated
      - TSS bam files
      - Full length bam files
  - Generate gene-level count data
      - From total bam file
      - From deduplicated bam file

The count data generated here will be mainly used for the QC purpose
only.

# Set up

``` r
processors <- 24

library("GenomicAlignments")
library("rtracklayer")
library("Rsubread")
library("kableExtra")

temp <- sapply(list.files("../functions", full.names = TRUE), source)
```

``` r
sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")
results.dir <- file.path("../../results")

## Input bam files
s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
star.aligned.bam.dir <- file.path(s2.alignment.dir, "s2-1-b-star-aligned_bam")

## Input annotation
annot.dir <- file.path("../../annotation/")
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")

all.tx.gtf <- file.path(
    annot.ps.dir,
    "all-transcript.gtf"
)

annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)


## Output file
s2.2.processed.bam.dir <-  file.path(s2.alignment.dir, "s2-2-processed-data")
s2.2.1.tss.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-1-tss-bam")
s2.2.2.dedup.tss.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-2-dedup-tss-bam")
s2.2.3.dedup.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-3-dedup-bam") 
s2.2.4.gene.count.dir <- file.path(s2.2.processed.bam.dir, "s2-2-4-gene-count")
s2.2.4.1.gene.count.total.dir <- file.path(s2.2.4.gene.count.dir, "s2-2-4-1-gene-count-total")
s2.2.4.2.gene.count.dedup.dir <- file.path(s2.2.4.gene.count.dir, "s2-2-4-2-gene-count-dedup")

create.dirs(c(
    s2.2.processed.bam.dir,
    s2.2.1.tss.bam.dir,
    s2.2.2.dedup.tss.bam.dir,
    s2.2.3.dedup.bam.dir,
    s2.2.4.gene.count.dir,
    s2.2.4.1.gene.count.total.dir,
    s2.2.4.2.gene.count.dedup.dir
))
```

# Generate TSS bam file

``` r
sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]

runGenerateTssBam <- function(sample.name, star.aligned.bam.dir, s2.2.1.tss.bam.dir){

    generateTssBam.rscript <- "./functions/generateTssBam.R"

    generateTssBam.cmd <- paste(
        "Rscript",
        generateTssBam.rscript,
        "-s", sample.name,
        "-b", star.aligned.bam.dir,
        "-o", s2.2.1.tss.bam.dir
    )

    system.cat(generateTssBam.cmd)
}


temp <- mclapply(
    sample.names,
    runGenerateTssBam,
    star.aligned.bam.dir = star.aligned.bam.dir,
    s2.2.1.tss.bam.dir = s2.2.1.tss.bam.dir,
    mc.cores = processors
)
```

# Deduplication with UMI

``` r
dedupTssBam <- function(sample.name, s2.2.1.tss.bam.dir, s2.2.2.dedup.tss.bam.dir){

    input.bam <- list.files(
        s2.2.1.tss.bam.dir,
        pattern = paste0(sample.name, ".tss.bam"),
        full.names = TRUE
    )

    output.bam <- file.path(
        s2.2.2.dedup.tss.bam.dir,
        gsub(".bam", ".dedup.tss.bam", basename(input.bam))
    )

    umi.dedup.cmd <- paste(
        "umi_tools",
        "dedup",
        "-I", input.bam,
        ## paste0("--output-stats=", dedup.stats.out),
        "-S", output.bam
    )

    dedup.out <- system.cat(umi.dedup.cmd)

    return(dedup.out)
}

dedup.outs <- mclapply(
    sample.names,
    dedupTssBam,
    s2.2.1.tss.bam.dir = s2.2.1.tss.bam.dir,
    s2.2.2.dedup.tss.bam.dir = s2.2.2.dedup.tss.bam.dir,
    mc.cores = processors
)
```

# Generate deduplicated bam file

``` r
runDedupBam <- function(sample.name, star.aligned.bam.dir, s2.2.2.dedup.tss.bam.dir, s2.2.3.dedup.bam.dir){

    dedupBam.rscript <- "./functions/dedupBam.R"

    dedupBam.cmd <- paste(
        "Rscript",
        dedupBam.rscript,
        "-s", sample.name,
        "-b", star.aligned.bam.dir,
        "-d", s2.2.2.dedup.tss.bam.dir,
        "-o", s2.2.3.dedup.bam.dir
    )

    system.cat(dedupBam.cmd)

    return()
}

temp <- mclapply(
    sample.names,
    runDedupBam,
    star.aligned.bam.dir = star.aligned.bam.dir,
    s2.2.2.dedup.tss.bam.dir = s2.2.2.dedup.tss.bam.dir,
    s2.2.3.dedup.bam.dir = s2.2.3.dedup.bam.dir,
    mc.cores = round(processors / 4)
)
```

# Create read count table by genes

``` r
primary.tx.essential.dt <- primary.tx.dt[, c("gene_id", "gene_name", "biotype"), with = FALSE]

primary.tx.essential.dt[, biotype := factor(
                              biotype, levels = c(
                                           "protein_coding",
                                           "lncRNA",
                                           "miRNA",
                                           "rRNA",
                                           "other"
                                       ))]

primary.tx.essential.dt <-
    primary.tx.essential.dt[order(biotype, decreasing = FALSE)][!duplicated(gene_id)]
setkey(primary.tx.essential.dt, gene_id)


geneCount <- function(sample.names, bam.dir, all.tx.gtf, processors){

    input.bams <- file.path(
        bam.dir,
        paste0(sample.names, ".bam")
    )

    fc <- featureCounts(
        input.bams,
        annot.ext = all.tx.gtf,
        isGTFAnnotationFile = TRUE,
        GTF.featureType = "exon",
        GTF.attrType = "gene_id",
        nthreads = processors,
        isPairedEnd = TRUE,
        strandSpecific = 1,
        primaryOnly = TRUE,
        requireBothEndsMapped = TRUE
    )

    count.dt <- data.table(fc$count, keep.rownames = TRUE)

    simplified.colnames <- gsub("\\/", "\\.", bam.dir) %>%
        gsub("\\-", "\\.", .) %>%
        gsub(., "", colnames(count.dt)) %>%
        gsub("^\\.", "", .) %>%
        gsub(".bam$", "", .) %>%
        gsub(".dedup$", "", .)
    
    setnames(count.dt, colnames(count.dt), simplified.colnames)
    setnames(count.dt, "rn", "gene_id")

    return(count.dt)
}

total.count.dt <- geneCount(
    sample.names,
    bam.dir = star.aligned.bam.dir,
    all.tx.gtf = all.tx.gtf,
    processors = processors
)
```

    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.2.1
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 206 BAM files                                    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0A.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0B.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo1.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo2.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo3.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo4.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo5.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo6.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo7.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo8.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0A.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0B.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo1.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo2.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo3.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo4.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo5.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo6.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo7.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo8.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0A.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0B.bam    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo1.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo2.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo3.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo4.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo5.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo6.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo7.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo8.bam     ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo1 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo2 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo3 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo4 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo5 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo6 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo7 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo8 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo1 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo2 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo3 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo4 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo5 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo6 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo7 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo8 ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0A.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0B.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo1.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo2.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo3.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo4.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo5.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo6.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo7.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo8.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0A.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0B.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo1.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo2.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo3.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo4.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo5.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo6.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo7.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo8.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0A.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0B.bam  ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo1.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo2.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo3.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo4.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo5.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo6.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo7.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo8.bam   ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo1.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo2.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo3.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo4.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo5.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo6.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo7.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo8.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo1.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo2.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo3.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo4.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo5.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo6.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo7.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo8.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo1.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo2.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo3.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo4.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo5.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo6.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo7.bam     ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo8.bam     ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo1.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo2.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo3.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo4.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo5.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo6.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo7.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo8.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo1.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo2.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo3.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo4.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo5.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo6.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo7.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo8.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo1.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo2.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo3.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo4.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo5.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo6.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo7.bam   ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo8.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo1.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo2.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo3.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo4.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo5.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo6.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo7.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo8.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo1.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo2.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo3.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo4.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo5.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo6.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo7.bam   ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo8.bam   ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo1 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo2 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo3 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo4 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo5 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo6 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo7 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo8 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo1 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo2 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo3 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo4 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo5 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo6 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo7 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo8 ... ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_1.bam                   ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_3.bam                   ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_4.bam                   ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_1.bam                 ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_3.bam                 ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_4.bam                 ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_1.bam                   ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_3.bam                   ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_4.bam                   ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_1.bam                 ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_3.bam                 ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_4.bam                 ||
    ## ||                           o total_786O_VHL_HIF1B_N_1.bam                   ||
    ## ||                           o total_786O_VHL_HIF1B_N_2.bam                   ||
    ## ||                           o total_786O_VHL_HIF1B_N_3.bam                   ||
    ## ||                           o total_786O_VHL_HIF1B_N_4.bam                   ||
    ## ||                           o total_786O_noVHL_HIF1B_N_1.bam                 ||
    ## ||                           o total_786O_noVHL_HIF1B_N_2.bam                 ||
    ## ||                           o total_786O_noVHL_HIF1B_N_3.bam                 ||
    ## ||                           o total_786O_noVHL_HIF1B_N_4.bam                 ||
    ## ||                           o total_786O_VHL_noHIF1B_N_1.bam                 ||
    ## ||                           o total_786O_VHL_noHIF1B_N_2.bam                 ||
    ## ||                           o total_786O_VHL_noHIF1B_N_3.bam                 ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_1.bam               ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_2.bam               ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_3.bam               ||
    ## ||                                                                            ||
    ## ||              Annotation : all-transcript.gtf (GTF)                         ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 24                                               ||
    ## ||                   Level : meta-feature level                               ||
    ## ||              Paired-end : yes                                              ||
    ## ||      Multimapping reads : counted                                          ||
    ## ||     Multiple alignments : primary alignment only                           ||
    ## || Multi-overlapping reads : not counted                                      ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## ||          Chimeric reads : counted                                          ||
    ## ||        Both ends mapped : required                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file all-transcript.gtf ...                                ||
    ## ||    Features : 1603934                                                      ||
    ## ||    Meta-features : 59617                                                   ||
    ## ||    Chromosomes/contigs : 25                                                ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0A.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 949164                                               ||
    ## ||    Successfully assigned alignments : 883975 (93.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0B.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4087485                                              ||
    ## ||    Successfully assigned alignments : 3871807 (94.7%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3173744                                              ||
    ## ||    Successfully assigned alignments : 3051158 (96.1%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2367573                                              ||
    ## ||    Successfully assigned alignments : 2282210 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3160002                                              ||
    ## ||    Successfully assigned alignments : 3049102 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4241775                                              ||
    ## ||    Successfully assigned alignments : 4115003 (97.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4983173                                              ||
    ## ||    Successfully assigned alignments : 4845749 (97.2%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3062270                                              ||
    ## ||    Successfully assigned alignments : 2983233 (97.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2203901                                              ||
    ## ||    Successfully assigned alignments : 2149070 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2727658                                              ||
    ## ||    Successfully assigned alignments : 2650525 (97.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0A.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1387740                                              ||
    ## ||    Successfully assigned alignments : 1303844 (94.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0B.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5055388                                              ||
    ## ||    Successfully assigned alignments : 4868647 (96.3%)                      ||
    ## ||    Running time : 0.11 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3401828                                              ||
    ## ||    Successfully assigned alignments : 3287530 (96.6%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2791233                                              ||
    ## ||    Successfully assigned alignments : 2684619 (96.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3444366                                              ||
    ## ||    Successfully assigned alignments : 3310991 (96.1%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3547100                                              ||
    ## ||    Successfully assigned alignments : 3443935 (97.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2182836                                              ||
    ## ||    Successfully assigned alignments : 2118791 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2232340                                              ||
    ## ||    Successfully assigned alignments : 2176024 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2565801                                              ||
    ## ||    Successfully assigned alignments : 2496061 (97.3%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2358297                                              ||
    ## ||    Successfully assigned alignments : 2297139 (97.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0A.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 999771                                               ||
    ## ||    Successfully assigned alignments : 926468 (92.7%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0B.bam...            ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5473635                                              ||
    ## ||    Successfully assigned alignments : 5198912 (95.0%)                      ||
    ## ||    Running time : 0.13 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4052946                                              ||
    ## ||    Successfully assigned alignments : 3872056 (95.5%)                      ||
    ## ||    Running time : 0.07 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2975321                                              ||
    ## ||    Successfully assigned alignments : 2856097 (96.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2646085                                              ||
    ## ||    Successfully assigned alignments : 2553341 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2755379                                              ||
    ## ||    Successfully assigned alignments : 2665104 (96.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2892451                                              ||
    ## ||    Successfully assigned alignments : 2817960 (97.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2458659                                              ||
    ## ||    Successfully assigned alignments : 2381840 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2244910                                              ||
    ## ||    Successfully assigned alignments : 2187165 (97.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1366080                                              ||
    ## ||    Successfully assigned alignments : 1335189 (97.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0A.bam...        ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 618458                                               ||
    ## ||    Successfully assigned alignments : 584437 (94.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0B.bam...        ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4825896                                              ||
    ## ||    Successfully assigned alignments : 4707473 (97.5%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo1.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4066319                                              ||
    ## ||    Successfully assigned alignments : 3963076 (97.5%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo2.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3783016                                              ||
    ## ||    Successfully assigned alignments : 3660593 (96.8%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo3.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2721873                                              ||
    ## ||    Successfully assigned alignments : 2612362 (96.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo4.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2009390                                              ||
    ## ||    Successfully assigned alignments : 1946268 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo5.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1018795                                              ||
    ## ||    Successfully assigned alignments : 981666 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo6.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 985247                                               ||
    ## ||    Successfully assigned alignments : 953089 (96.7%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo7.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 866040                                               ||
    ## ||    Successfully assigned alignments : 834289 (96.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo8.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1291527                                              ||
    ## ||    Successfully assigned alignments : 1246079 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0A.bam...        ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1134133                                              ||
    ## ||    Successfully assigned alignments : 1074989 (94.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0B.bam...        ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4077636                                              ||
    ## ||    Successfully assigned alignments : 3977680 (97.5%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo1.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7496076                                              ||
    ## ||    Successfully assigned alignments : 7284026 (97.2%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo2.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3398713                                              ||
    ## ||    Successfully assigned alignments : 3279947 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo3.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1506508                                              ||
    ## ||    Successfully assigned alignments : 1453859 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo4.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1226596                                              ||
    ## ||    Successfully assigned alignments : 1179928 (96.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo5.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1286464                                              ||
    ## ||    Successfully assigned alignments : 1245742 (96.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo6.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1116974                                              ||
    ## ||    Successfully assigned alignments : 1069981 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo7.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 861279                                               ||
    ## ||    Successfully assigned alignments : 831748 (96.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo8.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 457476                                               ||
    ## ||    Successfully assigned alignments : 442896 (96.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0A.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 861125                                               ||
    ## ||    Successfully assigned alignments : 796032 (92.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0B.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5902895                                              ||
    ## ||    Successfully assigned alignments : 5687404 (96.3%)                      ||
    ## ||    Running time : 0.11 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5556612                                              ||
    ## ||    Successfully assigned alignments : 5371820 (96.7%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4674751                                              ||
    ## ||    Successfully assigned alignments : 4514602 (96.6%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3245917                                              ||
    ## ||    Successfully assigned alignments : 3141300 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2275890                                              ||
    ## ||    Successfully assigned alignments : 2203821 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2919694                                              ||
    ## ||    Successfully assigned alignments : 2845970 (97.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2118807                                              ||
    ## ||    Successfully assigned alignments : 2055233 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2048599                                              ||
    ## ||    Successfully assigned alignments : 1997600 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1638796                                              ||
    ## ||    Successfully assigned alignments : 1601972 (97.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0A.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 987390                                               ||
    ## ||    Successfully assigned alignments : 926191 (93.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0B.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6620789                                              ||
    ## ||    Successfully assigned alignments : 6366028 (96.2%)                      ||
    ## ||    Running time : 0.15 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3697792                                              ||
    ## ||    Successfully assigned alignments : 3569202 (96.5%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3065670                                              ||
    ## ||    Successfully assigned alignments : 2954648 (96.4%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2837458                                              ||
    ## ||    Successfully assigned alignments : 2734277 (96.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3759840                                              ||
    ## ||    Successfully assigned alignments : 3632607 (96.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2724721                                              ||
    ## ||    Successfully assigned alignments : 2647790 (97.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2454163                                              ||
    ## ||    Successfully assigned alignments : 2387741 (97.3%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2105332                                              ||
    ## ||    Successfully assigned alignments : 2055438 (97.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1880570                                              ||
    ## ||    Successfully assigned alignments : 1833647 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0A.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2008054                                              ||
    ## ||    Successfully assigned alignments : 1879039 (93.6%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0B.bam...          ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4758340                                              ||
    ## ||    Successfully assigned alignments : 4640588 (97.5%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5549486                                              ||
    ## ||    Successfully assigned alignments : 5383132 (97.0%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2708370                                              ||
    ## ||    Successfully assigned alignments : 2625534 (96.9%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2737933                                              ||
    ## ||    Successfully assigned alignments : 2642361 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2841244                                              ||
    ## ||    Successfully assigned alignments : 2750780 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2868079                                              ||
    ## ||    Successfully assigned alignments : 2776074 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1780097                                              ||
    ## ||    Successfully assigned alignments : 1726271 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1539420                                              ||
    ## ||    Successfully assigned alignments : 1494524 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1058895                                              ||
    ## ||    Successfully assigned alignments : 1025175 (96.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0A.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1050209                                              ||
    ## ||    Successfully assigned alignments : 997335 (95.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0B.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4887915                                              ||
    ## ||    Successfully assigned alignments : 4774757 (97.7%)                      ||
    ## ||    Running time : 0.09 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo1.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6291519                                              ||
    ## ||    Successfully assigned alignments : 6116346 (97.2%)                      ||
    ## ||    Running time : 0.08 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo2.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3391530                                              ||
    ## ||    Successfully assigned alignments : 3280844 (96.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo3.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2590664                                              ||
    ## ||    Successfully assigned alignments : 2499700 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo4.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1952369                                              ||
    ## ||    Successfully assigned alignments : 1878686 (96.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo5.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1191003                                              ||
    ## ||    Successfully assigned alignments : 1149755 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo6.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 840573                                               ||
    ## ||    Successfully assigned alignments : 809358 (96.3%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo7.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 552841                                               ||
    ## ||    Successfully assigned alignments : 534636 (96.7%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo8.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 530711                                               ||
    ## ||    Successfully assigned alignments : 512737 (96.6%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0A.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 738086                                               ||
    ## ||    Successfully assigned alignments : 701824 (95.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0B.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4343297                                              ||
    ## ||    Successfully assigned alignments : 4216302 (97.1%)                      ||
    ## ||    Running time : 0.07 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo1.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5907972                                              ||
    ## ||    Successfully assigned alignments : 5771552 (97.7%)                      ||
    ## ||    Running time : 0.07 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo2.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2869773                                              ||
    ## ||    Successfully assigned alignments : 2792335 (97.3%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo3.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1574235                                              ||
    ## ||    Successfully assigned alignments : 1520597 (96.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo4.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1216467                                              ||
    ## ||    Successfully assigned alignments : 1172741 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo5.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 830606                                               ||
    ## ||    Successfully assigned alignments : 796207 (95.9%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo6.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 376767                                               ||
    ## ||    Successfully assigned alignments : 360164 (95.6%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo7.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 264474                                               ||
    ## ||    Successfully assigned alignments : 251244 (95.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo8.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 390027                                               ||
    ## ||    Successfully assigned alignments : 375413 (96.3%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2368860                                              ||
    ## ||    Successfully assigned alignments : 2262324 (95.5%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2129712                                              ||
    ## ||    Successfully assigned alignments : 2054852 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2748678                                              ||
    ## ||    Successfully assigned alignments : 2663462 (96.9%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2837201                                              ||
    ## ||    Successfully assigned alignments : 2762708 (97.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3667256                                              ||
    ## ||    Successfully assigned alignments : 3579305 (97.6%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2735000                                              ||
    ## ||    Successfully assigned alignments : 2668486 (97.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2632968                                              ||
    ## ||    Successfully assigned alignments : 2574722 (97.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2377541                                              ||
    ## ||    Successfully assigned alignments : 2324496 (97.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 520264                                               ||
    ## ||    Successfully assigned alignments : 496822 (95.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 880285                                               ||
    ## ||    Successfully assigned alignments : 850756 (96.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 980769                                               ||
    ## ||    Successfully assigned alignments : 949431 (96.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1546077                                              ||
    ## ||    Successfully assigned alignments : 1509831 (97.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 835485                                               ||
    ## ||    Successfully assigned alignments : 815269 (97.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 639620                                               ||
    ## ||    Successfully assigned alignments : 626903 (98.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 900577                                               ||
    ## ||    Successfully assigned alignments : 881023 (97.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1058100                                              ||
    ## ||    Successfully assigned alignments : 1035126 (97.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo1.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1363299                                              ||
    ## ||    Successfully assigned alignments : 1299840 (95.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo2.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 926287                                               ||
    ## ||    Successfully assigned alignments : 890012 (96.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo3.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 922874                                               ||
    ## ||    Successfully assigned alignments : 893296 (96.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo4.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 989007                                               ||
    ## ||    Successfully assigned alignments : 958892 (97.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo5.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1272322                                              ||
    ## ||    Successfully assigned alignments : 1242243 (97.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo6.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1021701                                              ||
    ## ||    Successfully assigned alignments : 993358 (97.2%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo7.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1034311                                              ||
    ## ||    Successfully assigned alignments : 1011206 (97.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo8.bam...             ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 668432                                               ||
    ## ||    Successfully assigned alignments : 654979 (98.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2448640                                              ||
    ## ||    Successfully assigned alignments : 2345297 (95.8%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2436589                                              ||
    ## ||    Successfully assigned alignments : 2357867 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2926763                                              ||
    ## ||    Successfully assigned alignments : 2839621 (97.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3019881                                              ||
    ## ||    Successfully assigned alignments : 2936117 (97.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3422363                                              ||
    ## ||    Successfully assigned alignments : 3339987 (97.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2635983                                              ||
    ## ||    Successfully assigned alignments : 2570870 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2353907                                              ||
    ## ||    Successfully assigned alignments : 2301171 (97.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2586346                                              ||
    ## ||    Successfully assigned alignments : 2529189 (97.8%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 748208                                               ||
    ## ||    Successfully assigned alignments : 723214 (96.7%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1002038                                              ||
    ## ||    Successfully assigned alignments : 965561 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1386294                                              ||
    ## ||    Successfully assigned alignments : 1345817 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1583698                                              ||
    ## ||    Successfully assigned alignments : 1541230 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1183396                                              ||
    ## ||    Successfully assigned alignments : 1156558 (97.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 714596                                               ||
    ## ||    Successfully assigned alignments : 698848 (97.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 625012                                               ||
    ## ||    Successfully assigned alignments : 612819 (98.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 741011                                               ||
    ## ||    Successfully assigned alignments : 726122 (98.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 705324                                               ||
    ## ||    Successfully assigned alignments : 672784 (95.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 542933                                               ||
    ## ||    Successfully assigned alignments : 524772 (96.7%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 889863                                               ||
    ## ||    Successfully assigned alignments : 862518 (96.9%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 886454                                               ||
    ## ||    Successfully assigned alignments : 862579 (97.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1285285                                              ||
    ## ||    Successfully assigned alignments : 1251441 (97.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 839310                                               ||
    ## ||    Successfully assigned alignments : 819025 (97.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 804056                                               ||
    ## ||    Successfully assigned alignments : 784163 (97.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 703317                                               ||
    ## ||    Successfully assigned alignments : 681544 (96.9%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1813597                                              ||
    ## ||    Successfully assigned alignments : 1737440 (95.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1779752                                              ||
    ## ||    Successfully assigned alignments : 1717942 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2107372                                              ||
    ## ||    Successfully assigned alignments : 2036634 (96.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2651042                                              ||
    ## ||    Successfully assigned alignments : 2579215 (97.3%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2424545                                              ||
    ## ||    Successfully assigned alignments : 2359201 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1552052                                              ||
    ## ||    Successfully assigned alignments : 1515938 (97.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1376314                                              ||
    ## ||    Successfully assigned alignments : 1342836 (97.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1477056                                              ||
    ## ||    Successfully assigned alignments : 1441673 (97.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo1.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3089435                                              ||
    ## ||    Successfully assigned alignments : 2936082 (95.0%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo2.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2901990                                              ||
    ## ||    Successfully assigned alignments : 2792104 (96.2%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo3.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2510251                                              ||
    ## ||    Successfully assigned alignments : 2429773 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo4.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3033210                                              ||
    ## ||    Successfully assigned alignments : 2942943 (97.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo5.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3338581                                              ||
    ## ||    Successfully assigned alignments : 3255948 (97.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo6.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2748408                                              ||
    ## ||    Successfully assigned alignments : 2671340 (97.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo7.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2685695                                              ||
    ## ||    Successfully assigned alignments : 2625581 (97.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo8.bam...           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1778794                                              ||
    ## ||    Successfully assigned alignments : 1740945 (97.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo1.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2167653                                              ||
    ## ||    Successfully assigned alignments : 2069072 (95.5%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo2.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2108517                                              ||
    ## ||    Successfully assigned alignments : 2034350 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo3.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2232515                                              ||
    ## ||    Successfully assigned alignments : 2167499 (97.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo4.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2447889                                              ||
    ## ||    Successfully assigned alignments : 2380376 (97.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo5.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2186887                                              ||
    ## ||    Successfully assigned alignments : 2138424 (97.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo6.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1369122                                              ||
    ## ||    Successfully assigned alignments : 1335318 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo7.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1392353                                              ||
    ## ||    Successfully assigned alignments : 1363440 (97.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo8.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1345210                                              ||
    ## ||    Successfully assigned alignments : 1317958 (98.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo1.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2923244                                              ||
    ## ||    Successfully assigned alignments : 2802025 (95.9%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo2.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2605508                                              ||
    ## ||    Successfully assigned alignments : 2523071 (96.8%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo3.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3792250                                              ||
    ## ||    Successfully assigned alignments : 3676438 (96.9%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo4.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4340189                                              ||
    ## ||    Successfully assigned alignments : 4227752 (97.4%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo5.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5208327                                              ||
    ## ||    Successfully assigned alignments : 5071859 (97.4%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo6.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3503248                                              ||
    ## ||    Successfully assigned alignments : 3419549 (97.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo7.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2645109                                              ||
    ## ||    Successfully assigned alignments : 2581075 (97.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo8.bam...         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2474250                                              ||
    ## ||    Successfully assigned alignments : 2411322 (97.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_1.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 13774334                                             ||
    ## ||    Successfully assigned alignments : 13131343 (95.3%)                     ||
    ## ||    Running time : 0.12 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_3.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 13796719                                             ||
    ## ||    Successfully assigned alignments : 13114357 (95.1%)                     ||
    ## ||    Running time : 0.15 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_4.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 13574123                                             ||
    ## ||    Successfully assigned alignments : 12926333 (95.2%)                     ||
    ## ||    Running time : 0.13 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_1.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 16130518                                             ||
    ## ||    Successfully assigned alignments : 15338752 (95.1%)                     ||
    ## ||    Running time : 0.19 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_3.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 12663798                                             ||
    ## ||    Successfully assigned alignments : 11949241 (94.4%)                     ||
    ## ||    Running time : 0.12 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_4.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 10483906                                             ||
    ## ||    Successfully assigned alignments : 9979229 (95.2%)                      ||
    ## ||    Running time : 0.11 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_1.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 18011169                                             ||
    ## ||    Successfully assigned alignments : 17304312 (96.1%)                     ||
    ## ||    Running time : 0.20 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_3.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 17579023                                             ||
    ## ||    Successfully assigned alignments : 16828696 (95.7%)                     ||
    ## ||    Running time : 0.19 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_4.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 12884150                                             ||
    ## ||    Successfully assigned alignments : 12373853 (96.0%)                     ||
    ## ||    Running time : 0.13 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_1.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 18182265                                             ||
    ## ||    Successfully assigned alignments : 17539487 (96.5%)                     ||
    ## ||    Running time : 0.22 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_3.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 12389080                                             ||
    ## ||    Successfully assigned alignments : 11936911 (96.4%)                     ||
    ## ||    Running time : 0.12 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_4.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 13489428                                             ||
    ## ||    Successfully assigned alignments : 12961282 (96.1%)                     ||
    ## ||    Running time : 0.15 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_1.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6929247                                              ||
    ## ||    Successfully assigned alignments : 6684265 (96.5%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_2.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7262290                                              ||
    ## ||    Successfully assigned alignments : 6917869 (95.3%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_3.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6591785                                              ||
    ## ||    Successfully assigned alignments : 6316801 (95.8%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_4.bam...                           ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 8219539                                              ||
    ## ||    Successfully assigned alignments : 7868516 (95.7%)                      ||
    ## ||    Running time : 0.07 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_1.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7067243                                              ||
    ## ||    Successfully assigned alignments : 6782655 (96.0%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_2.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6416086                                              ||
    ## ||    Successfully assigned alignments : 6139245 (95.7%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_3.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5683044                                              ||
    ## ||    Successfully assigned alignments : 5452815 (95.9%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_4.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5767667                                              ||
    ## ||    Successfully assigned alignments : 5531826 (95.9%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_1.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6701908                                              ||
    ## ||    Successfully assigned alignments : 6447827 (96.2%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_2.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6170926                                              ||
    ## ||    Successfully assigned alignments : 5943453 (96.3%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_3.bam...                         ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5307077                                              ||
    ## ||    Successfully assigned alignments : 5094094 (96.0%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_1.bam...                       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 10137850                                             ||
    ## ||    Successfully assigned alignments : 9714974 (95.8%)                      ||
    ## ||    Running time : 0.08 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_2.bam...                       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 8197438                                              ||
    ## ||    Successfully assigned alignments : 7845978 (95.7%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_3.bam...                       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 9372800                                              ||
    ## ||    Successfully assigned alignments : 8997492 (96.0%)                      ||
    ## ||    Running time : 0.08 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

``` r
total.count.dt <- merge(primary.tx.essential.dt, total.count.dt)
total.count.file <- file.path(s2.2.4.1.gene.count.total.dir, "total_gene_count_table.csv")

fwrite(total.count.dt, total.count.file)

## Caculate gene count for dedup data

dedup.count.dt <- geneCount(
    paste0(sample.names, ".dedup"),
    bam.dir = s2.2.3.dedup.bam.dir,
    all.tx.gtf = all.tx.gtf,
    processors = processors
)
```

    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.2.1
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 206 BAM files                                    ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0A.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0B.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo1.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo2.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo3.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo4.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo5.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo6.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo7.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo8.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0A.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0B.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo1.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo2.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo3.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo4.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo5.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo6.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo7.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo8.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0A.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0B.de ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo1.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo2.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo3.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo4.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo5.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo6.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo7.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo8.ded ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo1 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo2 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo3 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo4 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo5 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo6 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo7 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo8 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo1 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo2 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo3 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo4 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo5 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo6 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo7 ... ||
    ## ||                           o polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo8 ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0A. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0B. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo1.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo2.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo3.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo4.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo5.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo6.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo7.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo8.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0A. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0B. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo1.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo2.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo3.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo4.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo5.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo6.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo7.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo8.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0A. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0B. ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo1.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo2.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo3.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo4.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo5.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo6.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo7.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo8.d ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_rib ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo1.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo2.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo3.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo4.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo5.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo6.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo7.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_1_NA_ribo8.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo1.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo2.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo3.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo4.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo5.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo6.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo7.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_2_NA_ribo8.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo1.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo2.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo3.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo4.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo5.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo6.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo7.ded ... ||
    ## ||                           o polysome_786O_VHL_EIF4E2_NA_4_NA_ribo8.ded ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo1.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo2.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo3.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo4.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo5.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo6.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo7.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo8.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo1.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo2.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo3.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo4.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo5.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo6.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo7.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo8.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo1.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo2.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo3.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo4.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo5.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo6.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo7.d ... ||
    ## ||                           o polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo8.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo1.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo2.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo3.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo4.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo5.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo6.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo7.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo8.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo1.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo2.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo3.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo4.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo5.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo6.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo7.d ... ||
    ## ||                           o polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo8.d ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo1 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo2 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo3 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo4 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo5 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo6 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo7 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo8 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo1 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo2 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo3 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo4 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo5 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo6 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo7 ... ||
    ## ||                           o polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo8 ... ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_1.dedup.bam             ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_3.dedup.bam             ||
    ## ||                           o total_RCC4_VHL_HIF1B_N_4.dedup.bam             ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_1.dedup.bam           ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_3.dedup.bam           ||
    ## ||                           o total_RCC4_noVHL_HIF1B_N_4.dedup.bam           ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_1.dedup.bam             ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_3.dedup.bam             ||
    ## ||                           o total_RCC4_VHL_HIF1B_H_4.dedup.bam             ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_1.dedup.bam           ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_3.dedup.bam           ||
    ## ||                           o total_RCC4_noVHL_HIF1B_H_4.dedup.bam           ||
    ## ||                           o total_786O_VHL_HIF1B_N_1.dedup.bam             ||
    ## ||                           o total_786O_VHL_HIF1B_N_2.dedup.bam             ||
    ## ||                           o total_786O_VHL_HIF1B_N_3.dedup.bam             ||
    ## ||                           o total_786O_VHL_HIF1B_N_4.dedup.bam             ||
    ## ||                           o total_786O_noVHL_HIF1B_N_1.dedup.bam           ||
    ## ||                           o total_786O_noVHL_HIF1B_N_2.dedup.bam           ||
    ## ||                           o total_786O_noVHL_HIF1B_N_3.dedup.bam           ||
    ## ||                           o total_786O_noVHL_HIF1B_N_4.dedup.bam           ||
    ## ||                           o total_786O_VHL_noHIF1B_N_1.dedup.bam           ||
    ## ||                           o total_786O_VHL_noHIF1B_N_2.dedup.bam           ||
    ## ||                           o total_786O_VHL_noHIF1B_N_3.dedup.bam           ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_1.dedup.bam         ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_2.dedup.bam         ||
    ## ||                           o total_786O_noVHL_noHIF1B_N_3.dedup.bam         ||
    ## ||                                                                            ||
    ## ||              Annotation : all-transcript.gtf (GTF)                         ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 24                                               ||
    ## ||                   Level : meta-feature level                               ||
    ## ||              Paired-end : yes                                              ||
    ## ||      Multimapping reads : counted                                          ||
    ## ||     Multiple alignments : primary alignment only                           ||
    ## || Multi-overlapping reads : not counted                                      ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## ||          Chimeric reads : counted                                          ||
    ## ||        Both ends mapped : required                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file all-transcript.gtf ...                                ||
    ## ||    Features : 1603934                                                      ||
    ## ||    Meta-features : 59617                                                   ||
    ## ||    Chromosomes/contigs : 25                                                ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0A.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 397534                                               ||
    ## ||    Successfully assigned alignments : 362871 (91.3%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0B.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1502698                                              ||
    ## ||    Successfully assigned alignments : 1393301 (92.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1521460                                              ||
    ## ||    Successfully assigned alignments : 1436386 (94.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1431069                                              ||
    ## ||    Successfully assigned alignments : 1364163 (95.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1953653                                              ||
    ## ||    Successfully assigned alignments : 1868642 (95.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2613547                                              ||
    ## ||    Successfully assigned alignments : 2518716 (96.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2923540                                              ||
    ## ||    Successfully assigned alignments : 2822848 (96.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1953006                                              ||
    ## ||    Successfully assigned alignments : 1891543 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1465903                                              ||
    ## ||    Successfully assigned alignments : 1422536 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1747956                                              ||
    ## ||    Successfully assigned alignments : 1692549 (96.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0A.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 515169                                               ||
    ## ||    Successfully assigned alignments : 474288 (92.1%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0B.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1602869                                              ||
    ## ||    Successfully assigned alignments : 1503618 (93.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1494533                                              ||
    ## ||    Successfully assigned alignments : 1416361 (94.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1647974                                              ||
    ## ||    Successfully assigned alignments : 1568115 (95.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2150500                                              ||
    ## ||    Successfully assigned alignments : 2048304 (95.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2227298                                              ||
    ## ||    Successfully assigned alignments : 2149603 (96.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1430854                                              ||
    ## ||    Successfully assigned alignments : 1381566 (96.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1498311                                              ||
    ## ||    Successfully assigned alignments : 1453446 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1691992                                              ||
    ## ||    Successfully assigned alignments : 1637057 (96.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1527737                                              ||
    ## ||    Successfully assigned alignments : 1480506 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0A.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 379073                                               ||
    ## ||    Successfully assigned alignments : 344046 (90.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0B.dedup.bam...      ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1828877                                              ||
    ## ||    Successfully assigned alignments : 1696511 (92.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1731929                                              ||
    ## ||    Successfully assigned alignments : 1621568 (93.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1742608                                              ||
    ## ||    Successfully assigned alignments : 1652070 (94.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1712033                                              ||
    ## ||    Successfully assigned alignments : 1640138 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1756471                                              ||
    ## ||    Successfully assigned alignments : 1688037 (96.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1886534                                              ||
    ## ||    Successfully assigned alignments : 1828729 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1657295                                              ||
    ## ||    Successfully assigned alignments : 1596911 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1531936                                              ||
    ## ||    Successfully assigned alignments : 1486230 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 963888                                               ||
    ## ||    Successfully assigned alignments : 938646 (97.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0A.dedup.bam...  ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 354268                                               ||
    ## ||    Successfully assigned alignments : 335873 (94.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0B.dedup.bam...  ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1754157                                              ||
    ## ||    Successfully assigned alignments : 1686780 (96.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo1.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2006996                                              ||
    ## ||    Successfully assigned alignments : 1935478 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo2.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2328857                                              ||
    ## ||    Successfully assigned alignments : 2237476 (96.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo3.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1879085                                              ||
    ## ||    Successfully assigned alignments : 1793342 (95.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo4.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1466990                                              ||
    ## ||    Successfully assigned alignments : 1414892 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo5.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 775198                                               ||
    ## ||    Successfully assigned alignments : 744348 (96.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo6.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 759890                                               ||
    ## ||    Successfully assigned alignments : 732463 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo7.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 673910                                               ||
    ## ||    Successfully assigned alignments : 646762 (96.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo8.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 935601                                               ||
    ## ||    Successfully assigned alignments : 898971 (96.1%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0A.dedup.bam...  ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 625198                                               ||
    ## ||    Successfully assigned alignments : 594594 (95.1%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0B.dedup.bam...  ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1613066                                              ||
    ## ||    Successfully assigned alignments : 1554084 (96.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo1.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3425705                                              ||
    ## ||    Successfully assigned alignments : 3293089 (96.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo2.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2175390                                              ||
    ## ||    Successfully assigned alignments : 2084845 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo3.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1127048                                              ||
    ## ||    Successfully assigned alignments : 1083119 (96.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo4.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 939400                                               ||
    ## ||    Successfully assigned alignments : 900590 (95.9%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo5.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1002799                                              ||
    ## ||    Successfully assigned alignments : 967990 (96.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo6.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 880780                                               ||
    ## ||    Successfully assigned alignments : 840793 (95.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo7.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 689516                                               ||
    ## ||    Successfully assigned alignments : 663868 (96.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo8.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 376756                                               ||
    ## ||    Successfully assigned alignments : 363902 (96.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0A.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 360538                                               ||
    ## ||    Successfully assigned alignments : 327392 (90.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0B.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1838157                                              ||
    ## ||    Successfully assigned alignments : 1720499 (93.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2210858                                              ||
    ## ||    Successfully assigned alignments : 2087645 (94.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2482305                                              ||
    ## ||    Successfully assigned alignments : 2361336 (95.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1995192                                              ||
    ## ||    Successfully assigned alignments : 1914867 (96.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1470917                                              ||
    ## ||    Successfully assigned alignments : 1416145 (96.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1882758                                              ||
    ## ||    Successfully assigned alignments : 1826045 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1428711                                              ||
    ## ||    Successfully assigned alignments : 1378644 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1382396                                              ||
    ## ||    Successfully assigned alignments : 1341973 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1104890                                              ||
    ## ||    Successfully assigned alignments : 1074999 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0A.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 414433                                               ||
    ## ||    Successfully assigned alignments : 382046 (92.2%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0B.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2028630                                              ||
    ## ||    Successfully assigned alignments : 1895258 (93.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1519424                                              ||
    ## ||    Successfully assigned alignments : 1430506 (94.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1686744                                              ||
    ## ||    Successfully assigned alignments : 1601432 (94.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1733442                                              ||
    ## ||    Successfully assigned alignments : 1653804 (95.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2262466                                              ||
    ## ||    Successfully assigned alignments : 2168536 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1774723                                              ||
    ## ||    Successfully assigned alignments : 1714759 (96.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1592791                                              ||
    ## ||    Successfully assigned alignments : 1540163 (96.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1389563                                              ||
    ## ||    Successfully assigned alignments : 1349141 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1261438                                              ||
    ## ||    Successfully assigned alignments : 1224499 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0A.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 763755                                               ||
    ## ||    Successfully assigned alignments : 698867 (91.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0B.dedup.bam...    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1351720                                              ||
    ## ||    Successfully assigned alignments : 1279272 (94.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2052628                                              ||
    ## ||    Successfully assigned alignments : 1946189 (94.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1440970                                              ||
    ## ||    Successfully assigned alignments : 1375983 (95.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1632675                                              ||
    ## ||    Successfully assigned alignments : 1559245 (95.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1798589                                              ||
    ## ||    Successfully assigned alignments : 1729794 (96.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1781844                                              ||
    ## ||    Successfully assigned alignments : 1713215 (96.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1182870                                              ||
    ## ||    Successfully assigned alignments : 1140131 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1044326                                              ||
    ## ||    Successfully assigned alignments : 1007906 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 731467                                               ||
    ## ||    Successfully assigned alignments : 705378 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0A.dedup.b ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 574976                                               ||
    ## ||    Successfully assigned alignments : 546452 (95.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0B.dedup.b ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1936059                                              ||
    ## ||    Successfully assigned alignments : 1866032 (96.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo1.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2922450                                              ||
    ## ||    Successfully assigned alignments : 2807021 (96.1%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo2.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2038911                                              ||
    ## ||    Successfully assigned alignments : 1954578 (95.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo3.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1715495                                              ||
    ## ||    Successfully assigned alignments : 1642808 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo4.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1373506                                              ||
    ## ||    Successfully assigned alignments : 1314272 (95.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo5.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 897127                                               ||
    ## ||    Successfully assigned alignments : 862322 (96.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo6.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 640892                                               ||
    ## ||    Successfully assigned alignments : 614419 (95.9%)                       ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo7.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 433575                                               ||
    ## ||    Successfully assigned alignments : 417799 (96.4%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo8.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 422270                                               ||
    ## ||    Successfully assigned alignments : 406797 (96.3%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0A.dedup.b ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 441574                                               ||
    ## ||    Successfully assigned alignments : 420461 (95.2%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0B.dedup.b ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1751976                                              ||
    ## ||    Successfully assigned alignments : 1678988 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo1.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2775154                                              ||
    ## ||    Successfully assigned alignments : 2682389 (96.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo2.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1750624                                              ||
    ## ||    Successfully assigned alignments : 1689032 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo3.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1098726                                              ||
    ## ||    Successfully assigned alignments : 1054530 (96.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo4.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 913695                                               ||
    ## ||    Successfully assigned alignments : 876976 (96.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo5.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 637438                                               ||
    ## ||    Successfully assigned alignments : 608361 (95.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo6.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 308197                                               ||
    ## ||    Successfully assigned alignments : 293554 (95.2%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo7.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 220892                                               ||
    ## ||    Successfully assigned alignments : 209163 (94.7%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo8.dedup.ba ... ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 310548                                               ||
    ## ||    Successfully assigned alignments : 297639 (95.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1204212                                              ||
    ## ||    Successfully assigned alignments : 1138301 (94.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1281456                                              ||
    ## ||    Successfully assigned alignments : 1226573 (95.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1703481                                              ||
    ## ||    Successfully assigned alignments : 1639952 (96.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1769227                                              ||
    ## ||    Successfully assigned alignments : 1714859 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2205521                                              ||
    ## ||    Successfully assigned alignments : 2141431 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1746638                                              ||
    ## ||    Successfully assigned alignments : 1695664 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1712925                                              ||
    ## ||    Successfully assigned alignments : 1667270 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_1_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1554654                                              ||
    ## ||    Successfully assigned alignments : 1513631 (97.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 320717                                               ||
    ## ||    Successfully assigned alignments : 304669 (95.0%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 596472                                               ||
    ## ||    Successfully assigned alignments : 572986 (96.1%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 698754                                               ||
    ## ||    Successfully assigned alignments : 673558 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1042720                                              ||
    ## ||    Successfully assigned alignments : 1014875 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 595331                                               ||
    ## ||    Successfully assigned alignments : 579478 (97.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 475323                                               ||
    ## ||    Successfully assigned alignments : 464796 (97.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 651854                                               ||
    ## ||    Successfully assigned alignments : 635946 (97.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_2_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 742524                                               ||
    ## ||    Successfully assigned alignments : 724185 (97.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo1.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 701918                                               ||
    ## ||    Successfully assigned alignments : 661313 (94.2%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo2.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 607909                                               ||
    ## ||    Successfully assigned alignments : 580482 (95.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo3.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 655317                                               ||
    ## ||    Successfully assigned alignments : 631675 (96.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo4.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 691257                                               ||
    ## ||    Successfully assigned alignments : 667947 (96.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo5.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 882028                                               ||
    ## ||    Successfully assigned alignments : 858620 (97.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo6.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 729143                                               ||
    ## ||    Successfully assigned alignments : 706487 (96.9%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo7.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 741891                                               ||
    ## ||    Successfully assigned alignments : 723041 (97.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_EIF4E2_NA_4_NA_ribo8.dedup.bam...       ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 498890                                               ||
    ## ||    Successfully assigned alignments : 487644 (97.7%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1228816                                              ||
    ## ||    Successfully assigned alignments : 1162894 (94.6%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1427297                                              ||
    ## ||    Successfully assigned alignments : 1367681 (95.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1776595                                              ||
    ## ||    Successfully assigned alignments : 1711065 (96.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1860296                                              ||
    ## ||    Successfully assigned alignments : 1798266 (96.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2059090                                              ||
    ## ||    Successfully assigned alignments : 1998341 (97.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1670679                                              ||
    ## ||    Successfully assigned alignments : 1620552 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1544260                                              ||
    ## ||    Successfully assigned alignments : 1502551 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1697591                                              ||
    ## ||    Successfully assigned alignments : 1652375 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 463060                                               ||
    ## ||    Successfully assigned alignments : 442264 (95.5%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 659331                                               ||
    ## ||    Successfully assigned alignments : 630815 (95.7%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 922433                                               ||
    ## ||    Successfully assigned alignments : 889908 (96.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1035424                                              ||
    ## ||    Successfully assigned alignments : 1002734 (96.8%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 821572                                               ||
    ## ||    Successfully assigned alignments : 800093 (97.4%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 513388                                               ||
    ## ||    Successfully assigned alignments : 500448 (97.5%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 461642                                               ||
    ## ||    Successfully assigned alignments : 451278 (97.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 547925                                               ||
    ## ||    Successfully assigned alignments : 535396 (97.7%)                       ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 433585                                               ||
    ## ||    Successfully assigned alignments : 411254 (94.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 381950                                               ||
    ## ||    Successfully assigned alignments : 367164 (96.1%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 622699                                               ||
    ## ||    Successfully assigned alignments : 600922 (96.5%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 631936                                               ||
    ## ||    Successfully assigned alignments : 613485 (97.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 856788                                               ||
    ## ||    Successfully assigned alignments : 831146 (97.0%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 594590                                               ||
    ## ||    Successfully assigned alignments : 578439 (97.3%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 577515                                               ||
    ## ||    Successfully assigned alignments : 561329 (97.2%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 510059                                               ||
    ## ||    Successfully assigned alignments : 493579 (96.8%)                       ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 939834                                               ||
    ## ||    Successfully assigned alignments : 889573 (94.7%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1081230                                              ||
    ## ||    Successfully assigned alignments : 1034833 (95.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1343050                                              ||
    ## ||    Successfully assigned alignments : 1289529 (96.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1671906                                              ||
    ## ||    Successfully assigned alignments : 1619352 (96.9%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1525628                                              ||
    ## ||    Successfully assigned alignments : 1478209 (96.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1047069                                              ||
    ## ||    Successfully assigned alignments : 1018958 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 939465                                               ||
    ## ||    Successfully assigned alignments : 913012 (97.2%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 997569                                               ||
    ## ||    Successfully assigned alignments : 969748 (97.2%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo1.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1538756                                              ||
    ## ||    Successfully assigned alignments : 1446267 (94.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo2.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1709654                                              ||
    ## ||    Successfully assigned alignments : 1629548 (95.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo3.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1620332                                              ||
    ## ||    Successfully assigned alignments : 1558177 (96.2%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo4.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1917462                                              ||
    ## ||    Successfully assigned alignments : 1851115 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo5.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2119765                                              ||
    ## ||    Successfully assigned alignments : 2058651 (97.1%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo6.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1813735                                              ||
    ## ||    Successfully assigned alignments : 1753878 (96.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo7.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1789538                                              ||
    ## ||    Successfully assigned alignments : 1741835 (97.3%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo8.dedup.bam...     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1233861                                              ||
    ## ||    Successfully assigned alignments : 1203194 (97.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo1.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1099389                                              ||
    ## ||    Successfully assigned alignments : 1037450 (94.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo2.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1260871                                              ||
    ## ||    Successfully assigned alignments : 1206461 (95.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo3.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1415475                                              ||
    ## ||    Successfully assigned alignments : 1366050 (96.5%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo4.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1525079                                              ||
    ## ||    Successfully assigned alignments : 1475279 (96.7%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo5.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1403829                                              ||
    ## ||    Successfully assigned alignments : 1366738 (97.4%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo6.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 927397                                               ||
    ## ||    Successfully assigned alignments : 900650 (97.1%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo7.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 958629                                               ||
    ## ||    Successfully assigned alignments : 935381 (97.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo8.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 928583                                               ||
    ## ||    Successfully assigned alignments : 906576 (97.6%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo1.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1539662                                              ||
    ## ||    Successfully assigned alignments : 1462597 (95.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo2.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1547936                                              ||
    ## ||    Successfully assigned alignments : 1485009 (95.9%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo3.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2248041                                              ||
    ## ||    Successfully assigned alignments : 2161548 (96.2%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo4.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2604689                                              ||
    ## ||    Successfully assigned alignments : 2521812 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo5.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2993238                                              ||
    ## ||    Successfully assigned alignments : 2896400 (96.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo6.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2182961                                              ||
    ## ||    Successfully assigned alignments : 2118295 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo7.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1730271                                              ||
    ## ||    Successfully assigned alignments : 1679015 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo8.dedup.bam...   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1622004                                              ||
    ## ||    Successfully assigned alignments : 1573551 (97.0%)                      ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_1.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7925684                                              ||
    ## ||    Successfully assigned alignments : 7419381 (93.6%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_3.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7703628                                              ||
    ## ||    Successfully assigned alignments : 7170505 (93.1%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_N_4.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7665733                                              ||
    ## ||    Successfully assigned alignments : 7156646 (93.4%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_1.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 8988734                                              ||
    ## ||    Successfully assigned alignments : 8355874 (93.0%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_3.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7394551                                              ||
    ## ||    Successfully assigned alignments : 6817571 (92.2%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_N_4.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5954679                                              ||
    ## ||    Successfully assigned alignments : 5548449 (93.2%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_1.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 9515954                                              ||
    ## ||    Successfully assigned alignments : 8986846 (94.4%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_3.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 9000165                                              ||
    ## ||    Successfully assigned alignments : 8437764 (93.8%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_VHL_HIF1B_H_4.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7087913                                              ||
    ## ||    Successfully assigned alignments : 6696669 (94.5%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_1.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 9671223                                              ||
    ## ||    Successfully assigned alignments : 9182018 (94.9%)                      ||
    ## ||    Running time : 0.06 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_3.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 6867605                                              ||
    ## ||    Successfully assigned alignments : 6513640 (94.8%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_RCC4_noVHL_HIF1B_H_4.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 7104152                                              ||
    ## ||    Successfully assigned alignments : 6695838 (94.3%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_1.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4421325                                              ||
    ## ||    Successfully assigned alignments : 4216076 (95.4%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_2.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4681190                                              ||
    ## ||    Successfully assigned alignments : 4391411 (93.8%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_3.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4249547                                              ||
    ## ||    Successfully assigned alignments : 4017850 (94.5%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_HIF1B_N_4.dedup.bam...                     ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5037900                                              ||
    ## ||    Successfully assigned alignments : 4746685 (94.2%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_1.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4323640                                              ||
    ## ||    Successfully assigned alignments : 4088543 (94.6%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_2.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4138233                                              ||
    ## ||    Successfully assigned alignments : 3903010 (94.3%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_3.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3713274                                              ||
    ## ||    Successfully assigned alignments : 3517393 (94.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_HIF1B_N_4.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3714253                                              ||
    ## ||    Successfully assigned alignments : 3516685 (94.7%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_1.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 4174755                                              ||
    ## ||    Successfully assigned alignments : 3966005 (95.0%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_2.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3883485                                              ||
    ## ||    Successfully assigned alignments : 3694123 (95.1%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_VHL_noHIF1B_N_3.dedup.bam...                   ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3319075                                              ||
    ## ||    Successfully assigned alignments : 3145050 (94.8%)                      ||
    ## ||    Running time : 0.03 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_1.dedup.bam...                 ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5734981                                              ||
    ## ||    Successfully assigned alignments : 5407205 (94.3%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_2.dedup.bam...                 ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5049582                                              ||
    ## ||    Successfully assigned alignments : 4763984 (94.3%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file total_786O_noVHL_noHIF1B_N_3.dedup.bam...                 ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5315539                                              ||
    ## ||    Successfully assigned alignments : 5017697 (94.4%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

``` r
dedup.count.dt <- merge(primary.tx.essential.dt, dedup.count.dt)
dedup.count.file <- file.path(s2.2.4.2.gene.count.dedup.dir, "dedup_gene_count_table.csv")

fwrite(dedup.count.dt, dedup.count.file)
```

# Session Information

``` r
sessionInfo()
```

    ## R version 4.0.0 (2020-04-24)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/five_prime_seq_for_VHL_loss_v0.1.1/lib/libopenblasp-r0.3.9.so
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
    ##  [1] knitr_1.28                  stringr_1.4.0              
    ##  [3] magrittr_1.5                data.table_1.12.8          
    ##  [5] dplyr_1.0.0                 khroma_1.3.0               
    ##  [7] ggplot2_3.3.1               kableExtra_1.1.0           
    ##  [9] Rsubread_2.2.1              rtracklayer_1.48.0         
    ## [11] GenomicAlignments_1.24.0    Rsamtools_2.4.0            
    ## [13] Biostrings_2.56.0           XVector_0.28.0             
    ## [15] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
    ## [17] matrixStats_0.56.0          Biobase_2.48.0             
    ## [19] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
    ## [21] IRanges_2.22.1              S4Vectors_0.26.0           
    ## [23] BiocGenerics_0.34.0         rmarkdown_2.2              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6           lattice_0.20-41        digest_0.6.25         
    ##  [4] R6_2.4.1               evaluate_0.14          httr_1.4.1            
    ##  [7] pillar_1.4.4           zlibbioc_1.34.0        rlang_0.4.6           
    ## [10] rstudioapi_0.11        Matrix_1.2-18          webshot_0.5.2         
    ## [13] BiocParallel_1.22.0    readr_1.3.1            RCurl_1.98-1.2        
    ## [16] munsell_0.5.0          compiler_4.0.0         xfun_0.14             
    ## [19] pkgconfig_2.0.3        htmltools_0.4.0        tidyselect_1.1.0      
    ## [22] tibble_3.0.1           GenomeInfoDbData_1.2.3 XML_3.99-0.3          
    ## [25] viridisLite_0.3.0      crayon_1.3.4           withr_2.2.0           
    ## [28] bitops_1.0-6           grid_4.0.0             gtable_0.3.0          
    ## [31] lifecycle_0.2.0        scales_1.1.1           stringi_1.4.6         
    ## [34] xml2_1.3.2             ellipsis_0.3.1         vctrs_0.3.1           
    ## [37] generics_0.0.2         tools_4.0.0            glue_1.4.1            
    ## [40] purrr_0.3.4            hms_0.5.3              yaml_2.2.1            
    ## [43] colorspace_1.4-1       rvest_0.3.5
