---
title: "s3-1 Alignment statistics"
author: "Yoichiro Sugimoto"
date: "12 November, 2021"
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
   html_document:
     highlight: haddock
     toc: yes
     toc_depth: 2
     keep_md: yes
     fig_width: 7
     fig_height: 7
---


# Overview

The number of reads by the RNA source is analysed.


```r
## Bioconductor packages
library("GenomicAlignments")
library("Biostrings")
library("GenomicFeatures")
## knitr
library("kableExtra")

## Specify the number of CPUs to be used
processors <- 8

temp <- sapply(list.files("../functions", full.names = TRUE), source)
```



```r
sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

## Input annotation
annot.dir <- file.path("../../annotation/")

annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
annot.R.file <- list.files(
    annot.ps.dir,
    pattern = glob2rx("*primary_transcript_annotation*.rdata"),
    full.names = TRUE
)
load(annot.R.file)

## Input files
lib2fqfilename.file <- file.path("../../data/sample_data/20210306_libname_to_fq_filename.csv")

results.dir <- file.path("../../results")

processed.fq.dir <- file.path(results.dir, "s1-processed_fastq")
processed.fq.step0.dir <- file.path(processed.fq.dir, "s1-1-Step0")
processed.fq.step4.0.dir <- file.path(processed.fq.dir, "s1-1-Step4-0")
processed.fq.step4.dir <- file.path(processed.fq.dir, "s1-1-Step4")
processed.fq.step4.2.dir <- file.path(processed.fq.dir, "s1-1-Step4-2-rRNA-ERCC")

s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
star.aligned.bam.dir <- file.path(s2.alignment.dir, "s2-1-b-star-aligned_bam")
s2.2.processed.bam.dir <-  file.path(s2.alignment.dir, "s2-2-processed-data")
s2.2.1.tss.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-1-tss-bam")
s2.2.2.dedup.tss.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-2-dedup-tss-bam")
s2.2.3.dedup.bam.dir <- file.path(s2.2.processed.bam.dir, "s2-2-3-dedup-bam") 
s2.2.4.gene.count.dir <- file.path(s2.2.processed.bam.dir, "s2-2-4-gene-count")
s2.2.4.1.gene.count.total.dir <- file.path(s2.2.4.gene.count.dir, "s2-2-4-1-gene-count-total")
s2.2.4.2.gene.count.dedup.dir <- file.path(s2.2.4.gene.count.dir, "s2-2-4-2-gene-count-dedup")

s3.alignment.stats.dir <- file.path(results.dir, "s3-alignment-statistics")

create.dirs(
    c(s3.alignment.stats.dir)
)
```




```r
sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]

total.count.file <- file.path(s2.2.4.1.gene.count.total.dir, "total_gene_count_table.csv")
total.count.dt <- fread(total.count.file)
```


# Functions to extract read numbers



```r
countReadsFromFq <- function(sample.name, fq.dir, use.fq.filename = FALSE){
 
    if(use.fq.filename){
        input.fq <- file.path(fq.dir, sample.name)
    } else {
        input.fq <- list.files(
            fq.dir,
            pattern = paste0(
                "(", sample.name, ".*R1.fastq.gz|",
                sample.name, ".*fastq.1.gz)"
            ),
            full.names = TRUE
        )
    }
    
    rnum <- as.integer(
            system.cat(paste("zcat", input.fq, "| echo $((`wc -l`/4))")) 
    )

    return(rnum)
}


countReadsFromTssBam <- function(sample.name, bam.dir){

    input.bam <- list.files(
        bam.dir,
        pattern = glob2rx(paste0(sample.name, "*.bam")),
        full.names = TRUE
    )

    rnum <- countBam(input.bam)$records

    return(rnum)
}
```

# Alignment statistics by the RNA source

## rRNA and ERCC spike-in RNA aligned read counts


```r
ercc.fa.file <- file.path("../../data/ERCC_sequence/SRM2374_putative_T7_products_NoPolyA_v1.fasta.txt")

ercc.seq <- readDNAStringSet(ercc.fa.file, format = "fasta")
ercc.count.ref.dt <- data.table(
   RNA_source = c(names(ercc.seq), "rRNA_prerRNA", "others")
)
ercc.count.ref.dt <- ercc.count.ref.dt[order(RNA_source)][
    RNA_source != "080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev"
]

countERCC.rRNA <- function(sample.name, processed.fq.step4.2.dir){

    ## print(sample.name)
    
    rrna.ercc.bam.file <- file.path(
        processed.fq.step4.2.dir,
        paste0(sample.name, ".sorted.bam")
    )

    rrna.ercc.bam <- readGAlignmentPairs(
        rrna.ercc.bam.file,
        param = ScanBamParam(what = "mapq"),
        use.names = TRUE
    )

    rrna.ercc.unique.dt <- rrna.ercc.bam %>%
        as.data.frame %>%
    data.table(keep.rownames = "rname")

    ## sanity check
    if(nrow(rrna.ercc.unique.dt) != nrow(rrna.ercc.unique.dt[!duplicated(rname)])){
        stop("Duplicated reads!!")
    } else {"OK"}

    rrna.ids <- unique(rrna.ercc.unique.dt[, seqnames.first]) %>%
        {.[!grepl("^ERCC-", .)]}
    
    ## For ERCC counting, only considering reads uniquly aligned
    rrna.ercc.unique.dt[, `:=`(
        counting_flag = case_when(
            seqnames.first == seqnames.last &
            grepl("^ERCC-", seqnames.first) &
            strand.first == "+" & strand.last == "-" &
            mapq.first > 5 & mapq.last > 5 ~ as.character(seqnames.first),
            strand.first == "+" & strand.last == "-" &
            seqnames.first == seqnames.last &
            seqnames.first %in% rrna.ids &
            seqnames.last %in% rrna.ids ~ "rRNA_prerRNA",
            TRUE ~ "others"
        )
    )]
        
    rrna.ercc.count.dt <- rrna.ercc.unique.dt[, .N, by = counting_flag]

    setnames(
        rrna.ercc.count.dt,
        old = c("counting_flag", "N"),
        new = c("RNA_source", sample.name)
    )
    
    return(rrna.ercc.count.dt)
}

rrna.ercc.count.dts <- mclapply(
    sample.names,
    countERCC.rRNA,
    processed.fq.step4.2.dir = processed.fq.step4.2.dir,
    mc.cores = processors
)

rrna.ercc.count.dt <- Reduce(
    function(...) merge(..., all = TRUE, by = "RNA_source"),
    c(list(ercc.count.ref.dt), rrna.ercc.count.dts)
)

setnafill(
    rrna.ercc.count.dt,
    fill = 0,
    cols = sample.names
)

ercc.count.dt <- copy(rrna.ercc.count.dt)[grepl("^ERCC-", RNA_source)]
setnames(ercc.count.dt, old = "RNA_source", new = "ERCC_id")

ercc.count.file <- file.path(
    s3.alignment.stats.dir,
    "ERCC_count_per_sample.csv"
)

fwrite(
    ercc.count.dt, file = ercc.count.file
)

file.path(
    s3.alignment.stats.dir,
    "rRNA_ERCC_count_per_sample.csv"
) %>%
    {fwrite(rrna.ercc.count.dt, file = .)}

rrna.ercc.count.dt[, simplified_RNA_source := str_split_fixed(RNA_source, "-", n = 2)[, 1]]

sum.rrna.ercc.count.dt <- rrna.ercc.count.dt[
  , lapply(.SD, sum), .SDcols = sample.names, by = simplified_RNA_source
]
```

## RNA source analysis


```r
## Step 3: All adapter trimmed reads
all.rnums <- mclapply(
    sample.names,
    countReadsFromFq,
    fq.dir = processed.fq.step4.0.dir,
    mc.cores = processors
) %>% unlist
names(all.rnums) <- sample.names

## Step 4: All reads after rRNA removals
no.rRNA.rnums <- mclapply(
    sample.names,
    countReadsFromFq,
    fq.dir = processed.fq.step4.dir,
    mc.cores = processors
) %>% unlist
names(no.rRNA.rnums) <- sample.names

## All uniquly aligned reads
nodedup.rnums <- mclapply(
    sample.names,
    countReadsFromTssBam,
    bam.dir = s2.2.1.tss.bam.dir,
    mc.cores = processors
) %>% unlist
names(nodedup.rnums) <- sample.names

file.path(
    s3.alignment.stats.dir,
    "all-objects.Rdata"
) %>%
    {save(
         all.rnums, no.rRNA.rnums, nodedup.rnums, sum.rrna.ercc.count.dt,
         total.count.dt,
         file = .
     )}

createReadNumberDt <- function(all.rnums, no.rRNA.rnums, sum.rrna.ercc.count.dt, nodedup.rnums, count.dt, sample.dt, s3.alignment.stats.dir, out.prefix = "", dedup.flag = FALSE){

    aligned.rnums <- count.dt[, lapply(.SD, sum), .SDcols = sample.names] %>%
        unlist

    ercc.rnums <- sum.rrna.ercc.count.dt[
        simplified_RNA_source == "ERCC", sample.names, with = FALSE
    ] %>% unlist

    rrna.rnums <- sum.rrna.ercc.count.dt[
        simplified_RNA_source == "rRNA_prerRNA", sample.names, with = FALSE
    ] %>% unlist

    
    ## Sanity check
    if(
        all(names(aligned.rnums) == names(no.rRNA.rnums)) &
        all(names(no.rRNA.rnums) == names(ercc.rnums)) &
        all(names(ercc.rnums) == names(rrna.rnums))
    ){
        "OK"
    } else {stop()}

    not.aligned.rnums <- all.rnums - aligned.rnums - ercc.rnums - rrna.rnums
    
    annotated.rnums.dt <- count.dt[
      , lapply(.SD, sum), by = list(biotype), .SDcols = sample.names
    ] %>%
        melt(
            id.vars = "biotype",
            variable.name = "sample_name"
        ) %>%
        dcast(
            sample_name ~ biotype
        )

    annotated.rnums.dt <- annotated.rnums.dt[order(match(sample_name, sample.names))]

    if(all(names(aligned.rnums) == annotated.rnums.dt[, sample_name])){
        "OK"
    } else {stop()}
    
    if(dedup.flag){
        ## (TODO) Unaligned vs deduplication is not well defined here.
        stop("dedup.flag = TRUE is not currently supported")
    } else {
        annotated.rnums.dt[, `:=`(
            rRNA_prerRNA = rrna.rnums,
            ERCC = ercc.rnums,
            not_aligned = not.aligned.rnums
        )]
    }

    ## rRNA annotation from rRNA alignment and gtf annotation merged as this is confusing
    annotated.rnums.dt[, rRNA := rRNA + rRNA_prerRNA]
    annotated.rnums.dt[, rRNA_prerRNA := NULL]
    
    type.ordered <- c(
        "protein_coding", "lncRNA", "miRNA", "other",
        "duplicated", "rRNA", "ERCC", "not_aligned"
    )

    annotated.rnums.dt <- annotated.rnums.dt[, c(
        "sample_name",
        type.ordered[type.ordered %in% colnames(annotated.rnums.dt)]
    ), with = FALSE]

    for.export.annotated.rnums.dt <- copy(annotated.rnums.dt)[
      , all_read := Reduce("+", .SD),
        .SDcols = type.ordered[type.ordered %in% colnames(annotated.rnums.dt)]
    ]

    for.export.annotated.rnums.dt <- merge(
        sample.dt[, .(
            sample_name, experiment, cell, VHL, HIF1B, EIF4E2, gRNA_id, oxygen, clone, treatment, fraction, input_volume
        )],
        for.export.annotated.rnums.dt,
        by = "sample_name"
    )

    for.export.annotated.rnums.dt <- for.export.annotated.rnums.dt[
        order(match(sample_name, sample.names))
    ]

    for.export.annotated.rnums.dt[, `:=`(
        experiment = case_when(
            experiment == "polysome" ~ "HP5",
            experiment == "total" ~ "5' end-Seq"
        )
    )]
    
    file.path(
        s3.alignment.stats.dir, paste0(out.prefix, "alignment-statistics.csv")
    ) %>%
        {fwrite(for.export.annotated.rnums.dt, file = .)}
    
    m.annotated.rnums.dt <- melt(
        annotated.rnums.dt,
        id.vars = "sample_name",
        variable.name = "read_type",
        value = "read_count"
    )

    m.annotated.rnums.dt[, read_type := factor(read_type, levels = rev(type.ordered))]

    m.annotated.rnums.dt[, `:=`(
        experiment_type = str_split_fixed(sample_name, "_", n = 5)[, 1],
        cell = str_split_fixed(sample_name, "_", n = 5)[, 2],
        oxygen = case_when(
            grepl("^total_", sample_name) ~ str_split_fixed(sample_name, "_", n = 6)[, 5]
        ) %>%
            {factor(., levels = c("N", "H"))},
        VHL = str_split_fixed(sample_name, "_", n = 5)[, 3] %>%
            {factor(., levels = c("noVHL", "VHL"))},
        treatment = case_when(
            str_split_fixed(sample_name, "_", n = 5)[, 1] == "total" ~ NA_character_,
            TRUE ~ str_split_fixed(sample_name, "_", n = 9)[, 7]
        ),
        fraction = case_when(
            str_split_fixed(sample_name, "_", n = 5)[, 1] == "total" ~ NA_character_,
            TRUE ~ str_split_fixed(sample_name, "_", n = 9)[, 8]
        ),
        HIF1BorEIF4E2 = str_split_fixed(sample_name, "_", n = 5)[, 4],
        short_sample_name = str_split_fixed(sample_name, "_", n = 3)[, 3] %>%
            {gsub("_NA", "", .)} %>%
            {gsub("(_HIF1B|_noHIF1B|_EIF4E2|noEIF4E2)", "", .)}
    )]

    m.annotated.rnums.dt <- merge(
        sample.dt[, .(sample_name, library_ID)],
        m.annotated.rnums.dt,
        by = "sample_name"
    ) 

    ## Summarize HP5 for RCC-4 VHL stats to report
    hp5.r4vhl.m.annotated.rnums.dt <- m.annotated.rnums.dt[
        experiment_type == "polysome" &
        cell == "RCC4" &
        VHL == "VHL" &
        treatment == "NA"
    ]

    hp5.r4vhl.m.annotated.rnums.dt[
      , read_count_per_library := sum(read_count), by = sample_name
    ]

    hp5.r4vhl.m.annotated.rnums.dt[, `:=`(
        count_ratio = read_count / read_count_per_library
    )]

    print("Average number of reads for HP5 (RCC-4 VHL)")
    hp5.r4vhl.m.annotated.rnums.dt[
        read_type == "protein_coding",
        list(
            read_number =
                mean(read_count_per_library) %>%
                prettyNum(big.mark = ",",scientific = FALSE)
        )
    ] %>% print

    print("Average number of reads per fraction (RCC-4 VHL)")
    hp5.r4vhl.m.annotated.rnums.dt[
        read_type == "protein_coding",
        list(
            read_number = mean(read_count_per_library)
        ), by = fraction
    ] %T>%
        print %>%
        {.[,
        list(
            read_number =
                mean(read_number) %>%
                prettyNum(big.mark = ",",scientific = FALSE)
        )           
        ]} %>% print

    print("Proportion of reads from mRNA per fraction (RCC-4 VHL)")
    hp5.r4vhl.m.annotated.rnums.dt[
        read_type == "protein_coding",
        list(
            ratio_of_reads = mean(count_ratio)
        ), by = fraction
    ] %T>%
        print %>%
        {.[,
           list(
               read_number =
                   mean(ratio_of_reads) %>%
                   {sprintf("%1.0f%%", 100*.)}
           )           
           ]} %>% print

    print("Proportion of reads from ERCC per fraction (RCC-4 VHL)")
    hp5.r4vhl.m.annotated.rnums.dt[
        read_type == "ERCC",
        list(
            ratio_of_reads = mean(count_ratio)
        ), by = fraction
    ] %T>%
        print %>%
        {.[,
           list(
               read_number =
                   mean(ratio_of_reads) %>%
                   {sprintf("%1.2f%%", 100*.)}
           )           
           ]} %>% print

    read.type.colors <- c(
        "protein_coding" = "#2166AC",
        "lncRNA" = "skyblue1",
        "miRNA" = "paleturquoise3",
        "other" = "paleturquoise2",
        "not_annotated" = "paleturquoise1",
        "duplicated" = "papayawhip",
        "not_aligned" = "gray85",
        "rRNA" = "gray26",
        "ERCC" = "darkolivegreen4"        
    )
    
    g1 <- hp5.r4vhl.m.annotated.rnums.dt[
        ,
        list(
            read_number = mean(read_count)
        ), by = list(fraction, read_type)
    ] %>%
        ggplot(
            aes(
                x = str_extract(fraction, "[[:digit:]]") %>%
                    {case_when(
                         . == "1" ~ paste0(., " ribosome"),
                         . == "8" ~ paste0(., "+ ribosomes"),
                         TRUE ~ paste0(., " ribosomes")
                    )},
                y = read_number / 10^6,
                fill = read_type
            )
        ) +
        geom_bar(stat = "identity") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            aspect.ratio = 2,
            legend.position = "bottom"
        ) +
        xlab("Fraction") +
        ylab("Number of reads [10^6 reads]") +
        scale_fill_manual(
            values = read.type.colors
        ) +
        ggtitle("The number of reads by their identity\n(HP5 of RCC-4 VHL)")

    print(g1)
        
    return(for.export.annotated.rnums.dt)
}

total.for.export.annotated.rnums.dt <- createReadNumberDt(
    all.rnums = all.rnums,
    no.rRNA.rnums = no.rRNA.rnums,
    sum.rrna.ercc.count.dt = sum.rrna.ercc.count.dt,
    nodedup.rnums = nodedup.rnums,
    count.dt = total.count.dt,
    sample.dt = sample.dt,
    s3.alignment.stats.dir = s3.alignment.stats.dir,
    out.prefix = "",
    dedup.flag = FALSE
)
```

```
## [1] "Average number of reads for HP5 (RCC-4 VHL)"
##    read_number
## 1:   3,490,490
## [1] "Average number of reads per fraction (RCC-4 VHL)"
##     fraction read_number
##  1:   ribo0A     1643093
##  2:   ribo0B     6642110
##  3:    ribo1     4348908
##  4:    ribo2     3086370
##  5:    ribo3     3437579
##  6:    ribo4     3912157
##  7:    ribo5     3744135
##  8:    ribo6     2927857
##  9:    ribo7     2685902
## 10:    ribo8     2476788
##    read_number
## 1:   3,490,490
## [1] "Proportion of reads from mRNA per fraction (RCC-4 VHL)"
##     fraction ratio_of_reads
##  1:   ribo0A      0.5753763
##  2:   ribo0B      0.6471005
##  3:    ribo1      0.7372843
##  4:    ribo2      0.8101363
##  5:    ribo3      0.8403227
##  6:    ribo4      0.8531309
##  7:    ribo5      0.8539350
##  8:    ribo6      0.8469202
##  9:    ribo7      0.8395160
## 10:    ribo8      0.8307743
##    read_number
## 1:         78%
## [1] "Proportion of reads from ERCC per fraction (RCC-4 VHL)"
##     fraction ratio_of_reads
##  1:   ribo0A   0.0004798832
##  2:   ribo0B   0.0009919191
##  3:    ribo1   0.0006356504
##  4:    ribo2   0.0009172954
##  5:    ribo3   0.0008878090
##  6:    ribo4   0.0008853225
##  7:    ribo5   0.0012026413
##  8:    ribo6   0.0013062497
##  9:    ribo7   0.0016556181
## 10:    ribo8   0.0013815000
##    read_number
## 1:       0.10%
```

![](s3-1-alignment-statistics_files/figure-html/alignment statistics by the read origins total-1.png)<!-- -->



```r
total.for.export.annotated.rnums.dt %>%
    kable(., format.args = list(big.mark = ',')) %>%
    kable_styling(bootstrap_options = c("striped", "condensed"))
```

<table class="table table-striped table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> sample_name </th>
   <th style="text-align:left;"> experiment </th>
   <th style="text-align:left;"> cell </th>
   <th style="text-align:left;"> VHL </th>
   <th style="text-align:left;"> HIF1B </th>
   <th style="text-align:left;"> EIF4E2 </th>
   <th style="text-align:left;"> gRNA_id </th>
   <th style="text-align:left;"> oxygen </th>
   <th style="text-align:right;"> clone </th>
   <th style="text-align:left;"> treatment </th>
   <th style="text-align:left;"> fraction </th>
   <th style="text-align:right;"> input_volume </th>
   <th style="text-align:right;"> protein_coding </th>
   <th style="text-align:right;"> lncRNA </th>
   <th style="text-align:right;"> miRNA </th>
   <th style="text-align:right;"> other </th>
   <th style="text-align:right;"> rRNA </th>
   <th style="text-align:right;"> ERCC </th>
   <th style="text-align:right;"> not_aligned </th>
   <th style="text-align:right;"> all_read </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 807,386 </td>
   <td style="text-align:right;"> 24,944 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 51,607 </td>
   <td style="text-align:right;"> 287,638 </td>
   <td style="text-align:right;"> 615 </td>
   <td style="text-align:right;"> 233,869 </td>
   <td style="text-align:right;"> 1,406,094 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 3,563,920 </td>
   <td style="text-align:right;"> 128,602 </td>
   <td style="text-align:right;"> 114 </td>
   <td style="text-align:right;"> 179,148 </td>
   <td style="text-align:right;"> 929,982 </td>
   <td style="text-align:right;"> 4,999 </td>
   <td style="text-align:right;"> 750,289 </td>
   <td style="text-align:right;"> 5,557,054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,867,159 </td>
   <td style="text-align:right;"> 115,780 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 68,121 </td>
   <td style="text-align:right;"> 273,356 </td>
   <td style="text-align:right;"> 3,243 </td>
   <td style="text-align:right;"> 437,526 </td>
   <td style="text-align:right;"> 3,765,251 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,186,395 </td>
   <td style="text-align:right;"> 74,383 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 21,389 </td>
   <td style="text-align:right;"> 126,914 </td>
   <td style="text-align:right;"> 2,506 </td>
   <td style="text-align:right;"> 235,913 </td>
   <td style="text-align:right;"> 2,647,527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,965,600 </td>
   <td style="text-align:right;"> 71,285 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 12,174 </td>
   <td style="text-align:right;"> 137,907 </td>
   <td style="text-align:right;"> 3,121 </td>
   <td style="text-align:right;"> 272,034 </td>
   <td style="text-align:right;"> 3,462,154 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 4,042,502 </td>
   <td style="text-align:right;"> 62,811 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 9,647 </td>
   <td style="text-align:right;"> 186,415 </td>
   <td style="text-align:right;"> 4,076 </td>
   <td style="text-align:right;"> 337,510 </td>
   <td style="text-align:right;"> 4,642,997 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 4,782,190 </td>
   <td style="text-align:right;"> 55,114 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 8,377 </td>
   <td style="text-align:right;"> 222,393 </td>
   <td style="text-align:right;"> 6,700 </td>
   <td style="text-align:right;"> 381,797 </td>
   <td style="text-align:right;"> 5,456,633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,950,031 </td>
   <td style="text-align:right;"> 28,650 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 4,523 </td>
   <td style="text-align:right;"> 177,247 </td>
   <td style="text-align:right;"> 4,495 </td>
   <td style="text-align:right;"> 232,131 </td>
   <td style="text-align:right;"> 3,397,103 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,125,951 </td>
   <td style="text-align:right;"> 19,762 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3,341 </td>
   <td style="text-align:right;"> 149,614 </td>
   <td style="text-align:right;"> 3,907 </td>
   <td style="text-align:right;"> 172,610 </td>
   <td style="text-align:right;"> 2,475,200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,615,147 </td>
   <td style="text-align:right;"> 28,723 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 6,633 </td>
   <td style="text-align:right;"> 174,667 </td>
   <td style="text-align:right;"> 3,659 </td>
   <td style="text-align:right;"> 219,830 </td>
   <td style="text-align:right;"> 3,048,675 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,199,063 </td>
   <td style="text-align:right;"> 34,744 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 69,985 </td>
   <td style="text-align:right;"> 370,661 </td>
   <td style="text-align:right;"> 1,053 </td>
   <td style="text-align:right;"> 333,104 </td>
   <td style="text-align:right;"> 2,008,660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,521,362 </td>
   <td style="text-align:right;"> 150,345 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 196,845 </td>
   <td style="text-align:right;"> 1,211,004 </td>
   <td style="text-align:right;"> 6,841 </td>
   <td style="text-align:right;"> 764,929 </td>
   <td style="text-align:right;"> 6,851,413 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,095,860 </td>
   <td style="text-align:right;"> 113,379 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 78,229 </td>
   <td style="text-align:right;"> 462,367 </td>
   <td style="text-align:right;"> 2,415 </td>
   <td style="text-align:right;"> 467,726 </td>
   <td style="text-align:right;"> 4,220,023 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,575,136 </td>
   <td style="text-align:right;"> 81,209 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 28,227 </td>
   <td style="text-align:right;"> 188,408 </td>
   <td style="text-align:right;"> 3,499 </td>
   <td style="text-align:right;"> 324,095 </td>
   <td style="text-align:right;"> 3,200,612 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,221,614 </td>
   <td style="text-align:right;"> 75,433 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 13,886 </td>
   <td style="text-align:right;"> 205,871 </td>
   <td style="text-align:right;"> 3,817 </td>
   <td style="text-align:right;"> 343,314 </td>
   <td style="text-align:right;"> 3,863,987 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,381,216 </td>
   <td style="text-align:right;"> 53,935 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 8,743 </td>
   <td style="text-align:right;"> 253,975 </td>
   <td style="text-align:right;"> 4,221 </td>
   <td style="text-align:right;"> 277,678 </td>
   <td style="text-align:right;"> 3,979,804 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,089,189 </td>
   <td style="text-align:right;"> 24,457 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 5,124 </td>
   <td style="text-align:right;"> 183,619 </td>
   <td style="text-align:right;"> 3,534 </td>
   <td style="text-align:right;"> 175,488 </td>
   <td style="text-align:right;"> 2,481,428 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,150,627 </td>
   <td style="text-align:right;"> 21,142 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 4,243 </td>
   <td style="text-align:right;"> 210,381 </td>
   <td style="text-align:right;"> 3,680 </td>
   <td style="text-align:right;"> 176,499 </td>
   <td style="text-align:right;"> 2,566,583 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,468,912 </td>
   <td style="text-align:right;"> 22,353 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 4,774 </td>
   <td style="text-align:right;"> 266,351 </td>
   <td style="text-align:right;"> 5,796 </td>
   <td style="text-align:right;"> 209,722 </td>
   <td style="text-align:right;"> 2,977,929 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,267,119 </td>
   <td style="text-align:right;"> 23,272 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 6,737 </td>
   <td style="text-align:right;"> 269,761 </td>
   <td style="text-align:right;"> 4,523 </td>
   <td style="text-align:right;"> 191,982 </td>
   <td style="text-align:right;"> 2,763,404 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 840,527 </td>
   <td style="text-align:right;"> 24,272 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 61,627 </td>
   <td style="text-align:right;"> 325,928 </td>
   <td style="text-align:right;"> 724 </td>
   <td style="text-align:right;"> 261,412 </td>
   <td style="text-align:right;"> 1,514,525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,811,824 </td>
   <td style="text-align:right;"> 157,036 </td>
   <td style="text-align:right;"> 113 </td>
   <td style="text-align:right;"> 229,933 </td>
   <td style="text-align:right;"> 1,285,058 </td>
   <td style="text-align:right;"> 8,102 </td>
   <td style="text-align:right;"> 1,025,798 </td>
   <td style="text-align:right;"> 7,517,864 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,627,854 </td>
   <td style="text-align:right;"> 144,048 </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:right;"> 100,064 </td>
   <td style="text-align:right;"> 553,058 </td>
   <td style="text-align:right;"> 2,396 </td>
   <td style="text-align:right;"> 633,943 </td>
   <td style="text-align:right;"> 5,061,449 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,728,802 </td>
   <td style="text-align:right;"> 92,902 </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 34,336 </td>
   <td style="text-align:right;"> 221,712 </td>
   <td style="text-align:right;"> 2,429 </td>
   <td style="text-align:right;"> 330,736 </td>
   <td style="text-align:right;"> 3,410,971 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,480,780 </td>
   <td style="text-align:right;"> 60,898 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 11,622 </td>
   <td style="text-align:right;"> 184,964 </td>
   <td style="text-align:right;"> 2,312 </td>
   <td style="text-align:right;"> 245,986 </td>
   <td style="text-align:right;"> 2,986,597 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,612,782 </td>
   <td style="text-align:right;"> 43,470 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 8,824 </td>
   <td style="text-align:right;"> 213,998 </td>
   <td style="text-align:right;"> 2,234 </td>
   <td style="text-align:right;"> 232,335 </td>
   <td style="text-align:right;"> 3,113,670 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,778,699 </td>
   <td style="text-align:right;"> 32,954 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 6,289 </td>
   <td style="text-align:right;"> 256,750 </td>
   <td style="text-align:right;"> 3,149 </td>
   <td style="text-align:right;"> 216,486 </td>
   <td style="text-align:right;"> 3,294,345 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,352,997 </td>
   <td style="text-align:right;"> 24,249 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4,582 </td>
   <td style="text-align:right;"> 217,199 </td>
   <td style="text-align:right;"> 3,276 </td>
   <td style="text-align:right;"> 217,569 </td>
   <td style="text-align:right;"> 2,819,884 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,163,302 </td>
   <td style="text-align:right;"> 19,935 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 3,915 </td>
   <td style="text-align:right;"> 222,587 </td>
   <td style="text-align:right;"> 3,756 </td>
   <td style="text-align:right;"> 191,069 </td>
   <td style="text-align:right;"> 2,604,577 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,317,473 </td>
   <td style="text-align:right;"> 13,824 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 3,883 </td>
   <td style="text-align:right;"> 175,046 </td>
   <td style="text-align:right;"> 2,116 </td>
   <td style="text-align:right;"> 105,934 </td>
   <td style="text-align:right;"> 1,618,285 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 545,047 </td>
   <td style="text-align:right;"> 10,843 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 28,539 </td>
   <td style="text-align:right;"> 171,938 </td>
   <td style="text-align:right;"> 389 </td>
   <td style="text-align:right;"> 120,033 </td>
   <td style="text-align:right;"> 876,797 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,503,047 </td>
   <td style="text-align:right;"> 123,718 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 80,655 </td>
   <td style="text-align:right;"> 1,259,385 </td>
   <td style="text-align:right;"> 2,784 </td>
   <td style="text-align:right;"> 674,838 </td>
   <td style="text-align:right;"> 6,644,472 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,816,541 </td>
   <td style="text-align:right;"> 93,778 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 52,714 </td>
   <td style="text-align:right;"> 814,584 </td>
   <td style="text-align:right;"> 1,098 </td>
   <td style="text-align:right;"> 450,014 </td>
   <td style="text-align:right;"> 5,228,767 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,554,519 </td>
   <td style="text-align:right;"> 78,011 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 28,025 </td>
   <td style="text-align:right;"> 384,678 </td>
   <td style="text-align:right;"> 1,645 </td>
   <td style="text-align:right;"> 385,517 </td>
   <td style="text-align:right;"> 4,432,430 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,544,565 </td>
   <td style="text-align:right;"> 53,711 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 14,037 </td>
   <td style="text-align:right;"> 248,880 </td>
   <td style="text-align:right;"> 2,430 </td>
   <td style="text-align:right;"> 284,160 </td>
   <td style="text-align:right;"> 3,147,825 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,903,243 </td>
   <td style="text-align:right;"> 34,981 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 8,015 </td>
   <td style="text-align:right;"> 182,643 </td>
   <td style="text-align:right;"> 2,162 </td>
   <td style="text-align:right;"> 164,246 </td>
   <td style="text-align:right;"> 2,295,314 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 961,120 </td>
   <td style="text-align:right;"> 15,544 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 4,983 </td>
   <td style="text-align:right;"> 113,056 </td>
   <td style="text-align:right;"> 2,597 </td>
   <td style="text-align:right;"> 94,844 </td>
   <td style="text-align:right;"> 1,192,163 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 935,702 </td>
   <td style="text-align:right;"> 13,662 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 3,714 </td>
   <td style="text-align:right;"> 113,024 </td>
   <td style="text-align:right;"> 2,612 </td>
   <td style="text-align:right;"> 88,460 </td>
   <td style="text-align:right;"> 1,157,184 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 820,255 </td>
   <td style="text-align:right;"> 10,826 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 3,190 </td>
   <td style="text-align:right;"> 94,836 </td>
   <td style="text-align:right;"> 2,836 </td>
   <td style="text-align:right;"> 77,568 </td>
   <td style="text-align:right;"> 1,009,527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_3_Torin1_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,218,849 </td>
   <td style="text-align:right;"> 18,524 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 8,687 </td>
   <td style="text-align:right;"> 196,731 </td>
   <td style="text-align:right;"> 2,600 </td>
   <td style="text-align:right;"> 122,651 </td>
   <td style="text-align:right;"> 1,568,058 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,010,225 </td>
   <td style="text-align:right;"> 18,529 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 46,218 </td>
   <td style="text-align:right;"> 261,698 </td>
   <td style="text-align:right;"> 736 </td>
   <td style="text-align:right;"> 219,885 </td>
   <td style="text-align:right;"> 1,557,305 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,810,432 </td>
   <td style="text-align:right;"> 96,622 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 70,592 </td>
   <td style="text-align:right;"> 1,189,360 </td>
   <td style="text-align:right;"> 3,268 </td>
   <td style="text-align:right;"> 510,663 </td>
   <td style="text-align:right;"> 5,680,968 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 7,014,995 </td>
   <td style="text-align:right;"> 184,868 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 84,079 </td>
   <td style="text-align:right;"> 1,440,162 </td>
   <td style="text-align:right;"> 2,039 </td>
   <td style="text-align:right;"> 877,664 </td>
   <td style="text-align:right;"> 9,603,883 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,177,981 </td>
   <td style="text-align:right;"> 74,646 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 27,268 </td>
   <td style="text-align:right;"> 438,878 </td>
   <td style="text-align:right;"> 2,342 </td>
   <td style="text-align:right;"> 359,241 </td>
   <td style="text-align:right;"> 4,080,400 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,413,858 </td>
   <td style="text-align:right;"> 32,231 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 7,743 </td>
   <td style="text-align:right;"> 184,670 </td>
   <td style="text-align:right;"> 1,529 </td>
   <td style="text-align:right;"> 147,076 </td>
   <td style="text-align:right;"> 1,787,132 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,150,024 </td>
   <td style="text-align:right;"> 23,554 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 6,338 </td>
   <td style="text-align:right;"> 132,374 </td>
   <td style="text-align:right;"> 2,924 </td>
   <td style="text-align:right;"> 117,870 </td>
   <td style="text-align:right;"> 1,433,095 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,219,373 </td>
   <td style="text-align:right;"> 21,461 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 4,885 </td>
   <td style="text-align:right;"> 147,151 </td>
   <td style="text-align:right;"> 3,439 </td>
   <td style="text-align:right;"> 107,052 </td>
   <td style="text-align:right;"> 1,503,383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,049,485 </td>
   <td style="text-align:right;"> 16,659 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 3,827 </td>
   <td style="text-align:right;"> 113,750 </td>
   <td style="text-align:right;"> 4,206 </td>
   <td style="text-align:right;"> 114,548 </td>
   <td style="text-align:right;"> 1,302,485 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 817,760 </td>
   <td style="text-align:right;"> 11,312 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2,661 </td>
   <td style="text-align:right;"> 89,367 </td>
   <td style="text-align:right;"> 4,444 </td>
   <td style="text-align:right;"> 75,915 </td>
   <td style="text-align:right;"> 1,001,471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_VHL_EIF4E2_NA_4_Torin1_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 434,025 </td>
   <td style="text-align:right;"> 6,371 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2,497 </td>
   <td style="text-align:right;"> 74,276 </td>
   <td style="text-align:right;"> 1,744 </td>
   <td style="text-align:right;"> 40,187 </td>
   <td style="text-align:right;"> 559,103 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 726,094 </td>
   <td style="text-align:right;"> 21,373 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 48,534 </td>
   <td style="text-align:right;"> 276,025 </td>
   <td style="text-align:right;"> 1,281 </td>
   <td style="text-align:right;"> 214,671 </td>
   <td style="text-align:right;"> 1,288,008 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 5,296,741 </td>
   <td style="text-align:right;"> 182,452 </td>
   <td style="text-align:right;"> 101 </td>
   <td style="text-align:right;"> 208,107 </td>
   <td style="text-align:right;"> 1,417,431 </td>
   <td style="text-align:right;"> 6,024 </td>
   <td style="text-align:right;"> 906,852 </td>
   <td style="text-align:right;"> 8,017,708 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 5,072,376 </td>
   <td style="text-align:right;"> 188,214 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 111,158 </td>
   <td style="text-align:right;"> 575,617 </td>
   <td style="text-align:right;"> 4,619 </td>
   <td style="text-align:right;"> 688,957 </td>
   <td style="text-align:right;"> 6,641,011 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 4,330,891 </td>
   <td style="text-align:right;"> 136,930 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 46,726 </td>
   <td style="text-align:right;"> 269,998 </td>
   <td style="text-align:right;"> 5,168 </td>
   <td style="text-align:right;"> 462,511 </td>
   <td style="text-align:right;"> 5,252,277 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 3,051,587 </td>
   <td style="text-align:right;"> 74,808 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 14,872 </td>
   <td style="text-align:right;"> 179,943 </td>
   <td style="text-align:right;"> 3,711 </td>
   <td style="text-align:right;"> 290,591 </td>
   <td style="text-align:right;"> 3,615,545 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,157,790 </td>
   <td style="text-align:right;"> 38,297 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 7,712 </td>
   <td style="text-align:right;"> 123,451 </td>
   <td style="text-align:right;"> 3,945 </td>
   <td style="text-align:right;"> 188,910 </td>
   <td style="text-align:right;"> 2,520,126 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,803,533 </td>
   <td style="text-align:right;"> 36,159 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 6,250 </td>
   <td style="text-align:right;"> 163,871 </td>
   <td style="text-align:right;"> 4,694 </td>
   <td style="text-align:right;"> 218,970 </td>
   <td style="text-align:right;"> 3,233,502 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,028,463 </td>
   <td style="text-align:right;"> 22,345 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 4,406 </td>
   <td style="text-align:right;"> 137,796 </td>
   <td style="text-align:right;"> 4,709 </td>
   <td style="text-align:right;"> 180,998 </td>
   <td style="text-align:right;"> 2,378,735 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,974,488 </td>
   <td style="text-align:right;"> 19,285 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 3,811 </td>
   <td style="text-align:right;"> 133,473 </td>
   <td style="text-align:right;"> 5,494 </td>
   <td style="text-align:right;"> 173,080 </td>
   <td style="text-align:right;"> 2,309,647 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,579,724 </td>
   <td style="text-align:right;"> 17,525 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4,711 </td>
   <td style="text-align:right;"> 151,973 </td>
   <td style="text-align:right;"> 3,518 </td>
   <td style="text-align:right;"> 124,368 </td>
   <td style="text-align:right;"> 1,881,831 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 844,397 </td>
   <td style="text-align:right;"> 24,590 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 57,175 </td>
   <td style="text-align:right;"> 328,715 </td>
   <td style="text-align:right;"> 1,077 </td>
   <td style="text-align:right;"> 208,654 </td>
   <td style="text-align:right;"> 1,464,636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 5,941,159 </td>
   <td style="text-align:right;"> 210,829 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 213,934 </td>
   <td style="text-align:right;"> 1,458,798 </td>
   <td style="text-align:right;"> 6,166 </td>
   <td style="text-align:right;"> 1,140,909 </td>
   <td style="text-align:right;"> 8,971,895 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,372,581 </td>
   <td style="text-align:right;"> 122,498 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 74,074 </td>
   <td style="text-align:right;"> 428,774 </td>
   <td style="text-align:right;"> 2,389 </td>
   <td style="text-align:right;"> 498,127 </td>
   <td style="text-align:right;"> 4,498,491 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,834,704 </td>
   <td style="text-align:right;"> 90,285 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 29,620 </td>
   <td style="text-align:right;"> 216,253 </td>
   <td style="text-align:right;"> 2,627 </td>
   <td style="text-align:right;"> 323,709 </td>
   <td style="text-align:right;"> 3,497,235 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,655,807 </td>
   <td style="text-align:right;"> 65,607 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 12,828 </td>
   <td style="text-align:right;"> 214,979 </td>
   <td style="text-align:right;"> 3,288 </td>
   <td style="text-align:right;"> 265,650 </td>
   <td style="text-align:right;"> 3,218,194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,558,624 </td>
   <td style="text-align:right;"> 63,388 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 10,543 </td>
   <td style="text-align:right;"> 266,410 </td>
   <td style="text-align:right;"> 4,378 </td>
   <td style="text-align:right;"> 331,628 </td>
   <td style="text-align:right;"> 4,235,023 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,608,555 </td>
   <td style="text-align:right;"> 32,872 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 6,339 </td>
   <td style="text-align:right;"> 205,338 </td>
   <td style="text-align:right;"> 3,740 </td>
   <td style="text-align:right;"> 221,156 </td>
   <td style="text-align:right;"> 3,078,024 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,358,605 </td>
   <td style="text-align:right;"> 24,571 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 4,547 </td>
   <td style="text-align:right;"> 212,906 </td>
   <td style="text-align:right;"> 4,201 </td>
   <td style="text-align:right;"> 190,903 </td>
   <td style="text-align:right;"> 2,795,750 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,031,452 </td>
   <td style="text-align:right;"> 20,083 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 3,889 </td>
   <td style="text-align:right;"> 223,195 </td>
   <td style="text-align:right;"> 3,786 </td>
   <td style="text-align:right;"> 165,797 </td>
   <td style="text-align:right;"> 2,448,216 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,810,167 </td>
   <td style="text-align:right;"> 18,533 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 4,938 </td>
   <td style="text-align:right;"> 196,099 </td>
   <td style="text-align:right;"> 3,601 </td>
   <td style="text-align:right;"> 149,618 </td>
   <td style="text-align:right;"> 2,182,965 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,721,063 </td>
   <td style="text-align:right;"> 49,293 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 108,589 </td>
   <td style="text-align:right;"> 570,085 </td>
   <td style="text-align:right;"> 1,419 </td>
   <td style="text-align:right;"> 463,839 </td>
   <td style="text-align:right;"> 2,914,377 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,352,086 </td>
   <td style="text-align:right;"> 154,203 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 134,231 </td>
   <td style="text-align:right;"> 1,087,254 </td>
   <td style="text-align:right;"> 3,552 </td>
   <td style="text-align:right;"> 613,519 </td>
   <td style="text-align:right;"> 6,344,906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 5,082,210 </td>
   <td style="text-align:right;"> 181,845 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 119,009 </td>
   <td style="text-align:right;"> 736,511 </td>
   <td style="text-align:right;"> 2,808 </td>
   <td style="text-align:right;"> 785,714 </td>
   <td style="text-align:right;"> 6,908,159 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,513,338 </td>
   <td style="text-align:right;"> 81,158 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 30,998 </td>
   <td style="text-align:right;"> 231,861 </td>
   <td style="text-align:right;"> 1,951 </td>
   <td style="text-align:right;"> 269,355 </td>
   <td style="text-align:right;"> 3,128,699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,557,340 </td>
   <td style="text-align:right;"> 66,611 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 18,374 </td>
   <td style="text-align:right;"> 218,479 </td>
   <td style="text-align:right;"> 3,123 </td>
   <td style="text-align:right;"> 260,537 </td>
   <td style="text-align:right;"> 3,124,500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,688,820 </td>
   <td style="text-align:right;"> 49,935 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 11,990 </td>
   <td style="text-align:right;"> 234,038 </td>
   <td style="text-align:right;"> 3,293 </td>
   <td style="text-align:right;"> 254,127 </td>
   <td style="text-align:right;"> 3,242,237 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,728,178 </td>
   <td style="text-align:right;"> 38,149 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 9,708 </td>
   <td style="text-align:right;"> 268,298 </td>
   <td style="text-align:right;"> 4,025 </td>
   <td style="text-align:right;"> 256,379 </td>
   <td style="text-align:right;"> 3,304,776 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,701,010 </td>
   <td style="text-align:right;"> 19,713 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 5,526 </td>
   <td style="text-align:right;"> 217,907 </td>
   <td style="text-align:right;"> 3,068 </td>
   <td style="text-align:right;"> 157,828 </td>
   <td style="text-align:right;"> 2,105,070 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,474,277 </td>
   <td style="text-align:right;"> 15,562 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4,672 </td>
   <td style="text-align:right;"> 201,328 </td>
   <td style="text-align:right;"> 4,333 </td>
   <td style="text-align:right;"> 139,737 </td>
   <td style="text-align:right;"> 1,839,922 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,009,279 </td>
   <td style="text-align:right;"> 11,675 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4,215 </td>
   <td style="text-align:right;"> 137,946 </td>
   <td style="text-align:right;"> 2,350 </td>
   <td style="text-align:right;"> 99,646 </td>
   <td style="text-align:right;"> 1,265,116 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 945,696 </td>
   <td style="text-align:right;"> 15,483 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 36,134 </td>
   <td style="text-align:right;"> 221,675 </td>
   <td style="text-align:right;"> 615 </td>
   <td style="text-align:right;"> 167,254 </td>
   <td style="text-align:right;"> 1,386,878 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,590,653 </td>
   <td style="text-align:right;"> 112,273 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 71,793 </td>
   <td style="text-align:right;"> 1,213,715 </td>
   <td style="text-align:right;"> 1,876 </td>
   <td style="text-align:right;"> 621,413 </td>
   <td style="text-align:right;"> 6,611,759 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 5,922,548 </td>
   <td style="text-align:right;"> 129,045 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 64,688 </td>
   <td style="text-align:right;"> 965,969 </td>
   <td style="text-align:right;"> 1,801 </td>
   <td style="text-align:right;"> 705,592 </td>
   <td style="text-align:right;"> 7,789,703 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,190,982 </td>
   <td style="text-align:right;"> 65,399 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 24,430 </td>
   <td style="text-align:right;"> 406,526 </td>
   <td style="text-align:right;"> 1,440 </td>
   <td style="text-align:right;"> 334,555 </td>
   <td style="text-align:right;"> 4,023,364 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,435,823 </td>
   <td style="text-align:right;"> 50,719 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 13,118 </td>
   <td style="text-align:right;"> 279,268 </td>
   <td style="text-align:right;"> 1,665 </td>
   <td style="text-align:right;"> 238,530 </td>
   <td style="text-align:right;"> 3,019,161 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,833,606 </td>
   <td style="text-align:right;"> 35,933 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 9,121 </td>
   <td style="text-align:right;"> 192,913 </td>
   <td style="text-align:right;"> 2,353 </td>
   <td style="text-align:right;"> 184,847 </td>
   <td style="text-align:right;"> 2,258,798 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,125,470 </td>
   <td style="text-align:right;"> 19,558 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 4,718 </td>
   <td style="text-align:right;"> 127,634 </td>
   <td style="text-align:right;"> 2,204 </td>
   <td style="text-align:right;"> 105,334 </td>
   <td style="text-align:right;"> 1,384,927 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 793,464 </td>
   <td style="text-align:right;"> 12,598 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 3,283 </td>
   <td style="text-align:right;"> 89,890 </td>
   <td style="text-align:right;"> 2,127 </td>
   <td style="text-align:right;"> 74,271 </td>
   <td style="text-align:right;"> 975,646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 525,624 </td>
   <td style="text-align:right;"> 7,322 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 1,679 </td>
   <td style="text-align:right;"> 64,687 </td>
   <td style="text-align:right;"> 3,471 </td>
   <td style="text-align:right;"> 46,764 </td>
   <td style="text-align:right;"> 649,557 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_3_Torin1_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 502,776 </td>
   <td style="text-align:right;"> 7,386 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 2,565 </td>
   <td style="text-align:right;"> 70,241 </td>
   <td style="text-align:right;"> 1,461 </td>
   <td style="text-align:right;"> 46,795 </td>
   <td style="text-align:right;"> 631,234 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0A </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0A </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 660,644 </td>
   <td style="text-align:right;"> 10,529 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 30,634 </td>
   <td style="text-align:right;"> 190,432 </td>
   <td style="text-align:right;"> 835 </td>
   <td style="text-align:right;"> 117,901 </td>
   <td style="text-align:right;"> 1,010,989 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo0B </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo0B </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,049,134 </td>
   <td style="text-align:right;"> 91,177 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 75,947 </td>
   <td style="text-align:right;"> 1,202,540 </td>
   <td style="text-align:right;"> 3,452 </td>
   <td style="text-align:right;"> 689,459 </td>
   <td style="text-align:right;"> 6,111,752 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 5,576,415 </td>
   <td style="text-align:right;"> 129,011 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 66,073 </td>
   <td style="text-align:right;"> 1,051,453 </td>
   <td style="text-align:right;"> 1,985 </td>
   <td style="text-align:right;"> 642,459 </td>
   <td style="text-align:right;"> 7,467,444 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,713,295 </td>
   <td style="text-align:right;"> 57,528 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 21,482 </td>
   <td style="text-align:right;"> 420,561 </td>
   <td style="text-align:right;"> 1,508 </td>
   <td style="text-align:right;"> 274,084 </td>
   <td style="text-align:right;"> 3,488,487 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,479,261 </td>
   <td style="text-align:right;"> 31,097 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 10,213 </td>
   <td style="text-align:right;"> 175,843 </td>
   <td style="text-align:right;"> 2,191 </td>
   <td style="text-align:right;"> 150,612 </td>
   <td style="text-align:right;"> 1,849,241 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,142,954 </td>
   <td style="text-align:right;"> 23,083 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 6,681 </td>
   <td style="text-align:right;"> 135,920 </td>
   <td style="text-align:right;"> 2,117 </td>
   <td style="text-align:right;"> 119,096 </td>
   <td style="text-align:right;"> 1,429,873 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 776,994 </td>
   <td style="text-align:right;"> 14,720 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4,480 </td>
   <td style="text-align:right;"> 96,014 </td>
   <td style="text-align:right;"> 3,017 </td>
   <td style="text-align:right;"> 87,688 </td>
   <td style="text-align:right;"> 982,925 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 352,410 </td>
   <td style="text-align:right;"> 6,063 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1,687 </td>
   <td style="text-align:right;"> 38,441 </td>
   <td style="text-align:right;"> 2,985 </td>
   <td style="text-align:right;"> 39,704 </td>
   <td style="text-align:right;"> 441,294 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 245,723 </td>
   <td style="text-align:right;"> 4,159 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1,361 </td>
   <td style="text-align:right;"> 27,252 </td>
   <td style="text-align:right;"> 3,328 </td>
   <td style="text-align:right;"> 31,034 </td>
   <td style="text-align:right;"> 312,858 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_RCC4_noVHL_EIF4E2_NA_4_Torin1_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> Torin1 </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 366,978 </td>
   <td style="text-align:right;"> 6,114 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2,316 </td>
   <td style="text-align:right;"> 52,487 </td>
   <td style="text-align:right;"> 2,223 </td>
   <td style="text-align:right;"> 38,420 </td>
   <td style="text-align:right;"> 468,542 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,102,285 </td>
   <td style="text-align:right;"> 86,566 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 73,438 </td>
   <td style="text-align:right;"> 354,425 </td>
   <td style="text-align:right;"> 3,775 </td>
   <td style="text-align:right;"> 400,517 </td>
   <td style="text-align:right;"> 3,021,039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,963,755 </td>
   <td style="text-align:right;"> 62,483 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 28,589 </td>
   <td style="text-align:right;"> 183,606 </td>
   <td style="text-align:right;"> 3,983 </td>
   <td style="text-align:right;"> 254,562 </td>
   <td style="text-align:right;"> 2,497,002 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,593,612 </td>
   <td style="text-align:right;"> 56,895 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 12,914 </td>
   <td style="text-align:right;"> 179,334 </td>
   <td style="text-align:right;"> 4,035 </td>
   <td style="text-align:right;"> 259,544 </td>
   <td style="text-align:right;"> 3,106,371 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,712,168 </td>
   <td style="text-align:right;"> 42,432 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 8,084 </td>
   <td style="text-align:right;"> 199,608 </td>
   <td style="text-align:right;"> 4,312 </td>
   <td style="text-align:right;"> 234,579 </td>
   <td style="text-align:right;"> 3,201,204 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 3,527,796 </td>
   <td style="text-align:right;"> 44,327 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 7,149 </td>
   <td style="text-align:right;"> 246,719 </td>
   <td style="text-align:right;"> 6,681 </td>
   <td style="text-align:right;"> 292,507 </td>
   <td style="text-align:right;"> 4,125,206 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,634,310 </td>
   <td style="text-align:right;"> 29,403 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 4,754 </td>
   <td style="text-align:right;"> 228,547 </td>
   <td style="text-align:right;"> 5,474 </td>
   <td style="text-align:right;"> 237,531 </td>
   <td style="text-align:right;"> 3,140,036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,546,145 </td>
   <td style="text-align:right;"> 24,601 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 3,967 </td>
   <td style="text-align:right;"> 235,407 </td>
   <td style="text-align:right;"> 7,276 </td>
   <td style="text-align:right;"> 232,199 </td>
   <td style="text-align:right;"> 3,049,602 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,296,351 </td>
   <td style="text-align:right;"> 22,989 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 5,141 </td>
   <td style="text-align:right;"> 256,520 </td>
   <td style="text-align:right;"> 3,710 </td>
   <td style="text-align:right;"> 200,840 </td>
   <td style="text-align:right;"> 2,785,564 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 461,897 </td>
   <td style="text-align:right;"> 19,574 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 15,345 </td>
   <td style="text-align:right;"> 87,107 </td>
   <td style="text-align:right;"> 1,405 </td>
   <td style="text-align:right;"> 86,074 </td>
   <td style="text-align:right;"> 671,407 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 814,947 </td>
   <td style="text-align:right;"> 25,949 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9,854 </td>
   <td style="text-align:right;"> 50,525 </td>
   <td style="text-align:right;"> 2,153 </td>
   <td style="text-align:right;"> 92,940 </td>
   <td style="text-align:right;"> 996,374 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 925,150 </td>
   <td style="text-align:right;"> 20,605 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3,671 </td>
   <td style="text-align:right;"> 41,319 </td>
   <td style="text-align:right;"> 1,966 </td>
   <td style="text-align:right;"> 90,574 </td>
   <td style="text-align:right;"> 1,083,290 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,482,979 </td>
   <td style="text-align:right;"> 23,995 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 2,846 </td>
   <td style="text-align:right;"> 66,185 </td>
   <td style="text-align:right;"> 3,238 </td>
   <td style="text-align:right;"> 111,394 </td>
   <td style="text-align:right;"> 1,690,646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 803,071 </td>
   <td style="text-align:right;"> 10,437 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1,752 </td>
   <td style="text-align:right;"> 42,331 </td>
   <td style="text-align:right;"> 1,943 </td>
   <td style="text-align:right;"> 62,139 </td>
   <td style="text-align:right;"> 921,681 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 618,955 </td>
   <td style="text-align:right;"> 7,068 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 878 </td>
   <td style="text-align:right;"> 33,849 </td>
   <td style="text-align:right;"> 2,261 </td>
   <td style="text-align:right;"> 47,992 </td>
   <td style="text-align:right;"> 711,004 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 870,800 </td>
   <td style="text-align:right;"> 8,889 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1,326 </td>
   <td style="text-align:right;"> 51,860 </td>
   <td style="text-align:right;"> 3,097 </td>
   <td style="text-align:right;"> 68,615 </td>
   <td style="text-align:right;"> 1,004,593 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_2_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,021,869 </td>
   <td style="text-align:right;"> 10,980 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2,271 </td>
   <td style="text-align:right;"> 80,647 </td>
   <td style="text-align:right;"> 2,094 </td>
   <td style="text-align:right;"> 79,657 </td>
   <td style="text-align:right;"> 1,197,524 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,201,016 </td>
   <td style="text-align:right;"> 53,837 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 44,966 </td>
   <td style="text-align:right;"> 273,057 </td>
   <td style="text-align:right;"> 2,436 </td>
   <td style="text-align:right;"> 227,111 </td>
   <td style="text-align:right;"> 1,802,444 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 843,904 </td>
   <td style="text-align:right;"> 30,828 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 15,273 </td>
   <td style="text-align:right;"> 99,907 </td>
   <td style="text-align:right;"> 1,691 </td>
   <td style="text-align:right;"> 117,492 </td>
   <td style="text-align:right;"> 1,109,102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 864,265 </td>
   <td style="text-align:right;"> 23,107 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 5,916 </td>
   <td style="text-align:right;"> 67,393 </td>
   <td style="text-align:right;"> 2,123 </td>
   <td style="text-align:right;"> 92,846 </td>
   <td style="text-align:right;"> 1,055,658 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 935,617 </td>
   <td style="text-align:right;"> 18,579 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 4,681 </td>
   <td style="text-align:right;"> 69,632 </td>
   <td style="text-align:right;"> 2,041 </td>
   <td style="text-align:right;"> 89,435 </td>
   <td style="text-align:right;"> 1,119,999 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,220,609 </td>
   <td style="text-align:right;"> 18,258 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3,369 </td>
   <td style="text-align:right;"> 91,066 </td>
   <td style="text-align:right;"> 2,892 </td>
   <td style="text-align:right;"> 104,728 </td>
   <td style="text-align:right;"> 1,440,928 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 978,778 </td>
   <td style="text-align:right;"> 12,398 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 2,167 </td>
   <td style="text-align:right;"> 74,510 </td>
   <td style="text-align:right;"> 3,077 </td>
   <td style="text-align:right;"> 94,807 </td>
   <td style="text-align:right;"> 1,165,751 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 998,324 </td>
   <td style="text-align:right;"> 11,010 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1,865 </td>
   <td style="text-align:right;"> 82,369 </td>
   <td style="text-align:right;"> 3,495 </td>
   <td style="text-align:right;"> 96,404 </td>
   <td style="text-align:right;"> 1,193,472 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_EIF4E2_NA_4_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 645,649 </td>
   <td style="text-align:right;"> 7,237 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2,090 </td>
   <td style="text-align:right;"> 71,711 </td>
   <td style="text-align:right;"> 1,418 </td>
   <td style="text-align:right;"> 56,590 </td>
   <td style="text-align:right;"> 784,697 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,191,283 </td>
   <td style="text-align:right;"> 96,092 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 57,884 </td>
   <td style="text-align:right;"> 292,839 </td>
   <td style="text-align:right;"> 3,556 </td>
   <td style="text-align:right;"> 384,734 </td>
   <td style="text-align:right;"> 3,026,424 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,258,325 </td>
   <td style="text-align:right;"> 74,955 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 24,556 </td>
   <td style="text-align:right;"> 173,061 </td>
   <td style="text-align:right;"> 3,493 </td>
   <td style="text-align:right;"> 252,940 </td>
   <td style="text-align:right;"> 2,787,360 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,761,846 </td>
   <td style="text-align:right;"> 66,200 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 11,546 </td>
   <td style="text-align:right;"> 175,665 </td>
   <td style="text-align:right;"> 3,308 </td>
   <td style="text-align:right;"> 260,379 </td>
   <td style="text-align:right;"> 3,278,967 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,877,562 </td>
   <td style="text-align:right;"> 50,492 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 8,038 </td>
   <td style="text-align:right;"> 220,543 </td>
   <td style="text-align:right;"> 3,519 </td>
   <td style="text-align:right;"> 263,644 </td>
   <td style="text-align:right;"> 3,423,821 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 3,288,954 </td>
   <td style="text-align:right;"> 44,800 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 6,195 </td>
   <td style="text-align:right;"> 217,809 </td>
   <td style="text-align:right;"> 5,199 </td>
   <td style="text-align:right;"> 282,375 </td>
   <td style="text-align:right;"> 3,845,366 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,536,377 </td>
   <td style="text-align:right;"> 30,318 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 4,158 </td>
   <td style="text-align:right;"> 196,142 </td>
   <td style="text-align:right;"> 5,094 </td>
   <td style="text-align:right;"> 223,170 </td>
   <td style="text-align:right;"> 2,995,274 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,273,565 </td>
   <td style="text-align:right;"> 24,363 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3,226 </td>
   <td style="text-align:right;"> 189,835 </td>
   <td style="text-align:right;"> 5,348 </td>
   <td style="text-align:right;"> 204,788 </td>
   <td style="text-align:right;"> 2,701,137 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,495,151 </td>
   <td style="text-align:right;"> 28,874 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 5,141 </td>
   <td style="text-align:right;"> 276,685 </td>
   <td style="text-align:right;"> 3,165 </td>
   <td style="text-align:right;"> 219,577 </td>
   <td style="text-align:right;"> 3,028,615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 673,071 </td>
   <td style="text-align:right;"> 35,065 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 15,065 </td>
   <td style="text-align:right;"> 49,011 </td>
   <td style="text-align:right;"> 1,260 </td>
   <td style="text-align:right;"> 76,943 </td>
   <td style="text-align:right;"> 850,426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 922,974 </td>
   <td style="text-align:right;"> 30,993 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11,581 </td>
   <td style="text-align:right;"> 68,792 </td>
   <td style="text-align:right;"> 1,959 </td>
   <td style="text-align:right;"> 110,840 </td>
   <td style="text-align:right;"> 1,147,150 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,310,928 </td>
   <td style="text-align:right;"> 30,284 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 4,594 </td>
   <td style="text-align:right;"> 57,579 </td>
   <td style="text-align:right;"> 2,541 </td>
   <td style="text-align:right;"> 112,140 </td>
   <td style="text-align:right;"> 1,518,076 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,512,007 </td>
   <td style="text-align:right;"> 25,861 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 3,346 </td>
   <td style="text-align:right;"> 67,044 </td>
   <td style="text-align:right;"> 2,551 </td>
   <td style="text-align:right;"> 125,145 </td>
   <td style="text-align:right;"> 1,735,967 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,139,099 </td>
   <td style="text-align:right;"> 15,486 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 1,956 </td>
   <td style="text-align:right;"> 53,722 </td>
   <td style="text-align:right;"> 2,715 </td>
   <td style="text-align:right;"> 90,488 </td>
   <td style="text-align:right;"> 1,303,476 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 689,352 </td>
   <td style="text-align:right;"> 8,256 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,239 </td>
   <td style="text-align:right;"> 37,615 </td>
   <td style="text-align:right;"> 2,177 </td>
   <td style="text-align:right;"> 54,323 </td>
   <td style="text-align:right;"> 792,962 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 605,205 </td>
   <td style="text-align:right;"> 6,703 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 905 </td>
   <td style="text-align:right;"> 41,790 </td>
   <td style="text-align:right;"> 2,216 </td>
   <td style="text-align:right;"> 47,618 </td>
   <td style="text-align:right;"> 704,442 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_2_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 715,862 </td>
   <td style="text-align:right;"> 8,515 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1,744 </td>
   <td style="text-align:right;"> 48,405 </td>
   <td style="text-align:right;"> 1,870 </td>
   <td style="text-align:right;"> 54,246 </td>
   <td style="text-align:right;"> 830,643 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 626,176 </td>
   <td style="text-align:right;"> 25,769 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 20,820 </td>
   <td style="text-align:right;"> 90,836 </td>
   <td style="text-align:right;"> 1,136 </td>
   <td style="text-align:right;"> 122,749 </td>
   <td style="text-align:right;"> 887,505 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 502,834 </td>
   <td style="text-align:right;"> 15,286 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6,645 </td>
   <td style="text-align:right;"> 49,155 </td>
   <td style="text-align:right;"> 781 </td>
   <td style="text-align:right;"> 58,629 </td>
   <td style="text-align:right;"> 633,335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 840,571 </td>
   <td style="text-align:right;"> 17,835 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4,106 </td>
   <td style="text-align:right;"> 60,812 </td>
   <td style="text-align:right;"> 1,366 </td>
   <td style="text-align:right;"> 80,112 </td>
   <td style="text-align:right;"> 1,004,808 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 846,676 </td>
   <td style="text-align:right;"> 13,312 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2,582 </td>
   <td style="text-align:right;"> 64,287 </td>
   <td style="text-align:right;"> 1,319 </td>
   <td style="text-align:right;"> 77,801 </td>
   <td style="text-align:right;"> 1,005,984 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,233,324 </td>
   <td style="text-align:right;"> 15,580 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2,531 </td>
   <td style="text-align:right;"> 97,291 </td>
   <td style="text-align:right;"> 2,033 </td>
   <td style="text-align:right;"> 113,209 </td>
   <td style="text-align:right;"> 1,463,973 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 808,756 </td>
   <td style="text-align:right;"> 8,977 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1,288 </td>
   <td style="text-align:right;"> 72,771 </td>
   <td style="text-align:right;"> 1,616 </td>
   <td style="text-align:right;"> 72,747 </td>
   <td style="text-align:right;"> 966,158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 774,988 </td>
   <td style="text-align:right;"> 7,848 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1,320 </td>
   <td style="text-align:right;"> 80,777 </td>
   <td style="text-align:right;"> 1,817 </td>
   <td style="text-align:right;"> 75,206 </td>
   <td style="text-align:right;"> 941,962 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_EIF4E2_NA_4_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 672,907 </td>
   <td style="text-align:right;"> 7,096 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1,536 </td>
   <td style="text-align:right;"> 97,023 </td>
   <td style="text-align:right;"> 1,001 </td>
   <td style="text-align:right;"> 72,448 </td>
   <td style="text-align:right;"> 852,016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,615,749 </td>
   <td style="text-align:right;"> 68,201 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 53,457 </td>
   <td style="text-align:right;"> 281,084 </td>
   <td style="text-align:right;"> 3,328 </td>
   <td style="text-align:right;"> 281,032 </td>
   <td style="text-align:right;"> 2,302,883 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,637,650 </td>
   <td style="text-align:right;"> 55,935 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 24,339 </td>
   <td style="text-align:right;"> 141,427 </td>
   <td style="text-align:right;"> 2,983 </td>
   <td style="text-align:right;"> 197,800 </td>
   <td style="text-align:right;"> 2,060,148 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,976,581 </td>
   <td style="text-align:right;"> 48,836 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 11,184 </td>
   <td style="text-align:right;"> 133,317 </td>
   <td style="text-align:right;"> 3,386 </td>
   <td style="text-align:right;"> 196,766 </td>
   <td style="text-align:right;"> 2,370,102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,526,316 </td>
   <td style="text-align:right;"> 45,551 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 7,324 </td>
   <td style="text-align:right;"> 170,484 </td>
   <td style="text-align:right;"> 4,011 </td>
   <td style="text-align:right;"> 207,874 </td>
   <td style="text-align:right;"> 2,961,579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,322,100 </td>
   <td style="text-align:right;"> 31,898 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 5,171 </td>
   <td style="text-align:right;"> 163,216 </td>
   <td style="text-align:right;"> 4,617 </td>
   <td style="text-align:right;"> 189,542 </td>
   <td style="text-align:right;"> 2,716,571 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,495,452 </td>
   <td style="text-align:right;"> 17,822 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 2,649 </td>
   <td style="text-align:right;"> 128,633 </td>
   <td style="text-align:right;"> 3,391 </td>
   <td style="text-align:right;"> 121,103 </td>
   <td style="text-align:right;"> 1,769,063 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,326,381 </td>
   <td style="text-align:right;"> 13,862 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2,585 </td>
   <td style="text-align:right;"> 126,126 </td>
   <td style="text-align:right;"> 4,549 </td>
   <td style="text-align:right;"> 112,877 </td>
   <td style="text-align:right;"> 1,586,387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g1_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,421,829 </td>
   <td style="text-align:right;"> 15,785 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 4,046 </td>
   <td style="text-align:right;"> 160,388 </td>
   <td style="text-align:right;"> 2,658 </td>
   <td style="text-align:right;"> 121,857 </td>
   <td style="text-align:right;"> 1,726,574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,716,742 </td>
   <td style="text-align:right;"> 113,972 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 105,311 </td>
   <td style="text-align:right;"> 682,992 </td>
   <td style="text-align:right;"> 3,603 </td>
   <td style="text-align:right;"> 534,541 </td>
   <td style="text-align:right;"> 4,157,217 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,657,000 </td>
   <td style="text-align:right;"> 85,797 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 49,268 </td>
   <td style="text-align:right;"> 338,092 </td>
   <td style="text-align:right;"> 3,939 </td>
   <td style="text-align:right;"> 337,658 </td>
   <td style="text-align:right;"> 3,471,790 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,357,536 </td>
   <td style="text-align:right;"> 55,825 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 16,394 </td>
   <td style="text-align:right;"> 214,016 </td>
   <td style="text-align:right;"> 2,985 </td>
   <td style="text-align:right;"> 227,548 </td>
   <td style="text-align:right;"> 2,874,321 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,884,032 </td>
   <td style="text-align:right;"> 48,544 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 10,343 </td>
   <td style="text-align:right;"> 265,141 </td>
   <td style="text-align:right;"> 4,633 </td>
   <td style="text-align:right;"> 255,423 </td>
   <td style="text-align:right;"> 3,468,139 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,203,326 </td>
   <td style="text-align:right;"> 43,244 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 9,359 </td>
   <td style="text-align:right;"> 319,945 </td>
   <td style="text-align:right;"> 4,849 </td>
   <td style="text-align:right;"> 254,430 </td>
   <td style="text-align:right;"> 3,835,170 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,634,022 </td>
   <td style="text-align:right;"> 31,198 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 6,106 </td>
   <td style="text-align:right;"> 286,904 </td>
   <td style="text-align:right;"> 5,727 </td>
   <td style="text-align:right;"> 229,986 </td>
   <td style="text-align:right;"> 3,193,957 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,594,172 </td>
   <td style="text-align:right;"> 26,292 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 5,101 </td>
   <td style="text-align:right;"> 292,730 </td>
   <td style="text-align:right;"> 6,504 </td>
   <td style="text-align:right;"> 210,208 </td>
   <td style="text-align:right;"> 3,135,022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_VHL_noEIF4E2_g2_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 1,716,402 </td>
   <td style="text-align:right;"> 18,322 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 6,210 </td>
   <td style="text-align:right;"> 271,654 </td>
   <td style="text-align:right;"> 2,932 </td>
   <td style="text-align:right;"> 139,252 </td>
   <td style="text-align:right;"> 2,154,782 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,932,448 </td>
   <td style="text-align:right;"> 82,149 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 54,450 </td>
   <td style="text-align:right;"> 295,646 </td>
   <td style="text-align:right;"> 3,267 </td>
   <td style="text-align:right;"> 330,698 </td>
   <td style="text-align:right;"> 2,698,682 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,950,342 </td>
   <td style="text-align:right;"> 61,120 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 22,869 </td>
   <td style="text-align:right;"> 153,704 </td>
   <td style="text-align:right;"> 3,250 </td>
   <td style="text-align:right;"> 224,270 </td>
   <td style="text-align:right;"> 2,415,571 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,112,713 </td>
   <td style="text-align:right;"> 47,209 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 7,558 </td>
   <td style="text-align:right;"> 128,195 </td>
   <td style="text-align:right;"> 2,837 </td>
   <td style="text-align:right;"> 186,267 </td>
   <td style="text-align:right;"> 2,484,796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,336,549 </td>
   <td style="text-align:right;"> 38,279 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 5,528 </td>
   <td style="text-align:right;"> 150,533 </td>
   <td style="text-align:right;"> 3,584 </td>
   <td style="text-align:right;"> 197,561 </td>
   <td style="text-align:right;"> 2,732,049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 2,108,012 </td>
   <td style="text-align:right;"> 26,844 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3,552 </td>
   <td style="text-align:right;"> 155,492 </td>
   <td style="text-align:right;"> 3,539 </td>
   <td style="text-align:right;"> 166,070 </td>
   <td style="text-align:right;"> 2,463,524 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,318,186 </td>
   <td style="text-align:right;"> 15,113 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 2,000 </td>
   <td style="text-align:right;"> 104,531 </td>
   <td style="text-align:right;"> 3,254 </td>
   <td style="text-align:right;"> 110,650 </td>
   <td style="text-align:right;"> 1,553,745 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,348,063 </td>
   <td style="text-align:right;"> 13,649 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1,718 </td>
   <td style="text-align:right;"> 116,530 </td>
   <td style="text-align:right;"> 3,761 </td>
   <td style="text-align:right;"> 113,709 </td>
   <td style="text-align:right;"> 1,597,436 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g1_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g1 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1,301,528 </td>
   <td style="text-align:right;"> 13,875 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 2,545 </td>
   <td style="text-align:right;"> 141,105 </td>
   <td style="text-align:right;"> 2,494 </td>
   <td style="text-align:right;"> 101,542 </td>
   <td style="text-align:right;"> 1,563,099 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo1 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo1 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,610,561 </td>
   <td style="text-align:right;"> 117,063 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 74,351 </td>
   <td style="text-align:right;"> 386,005 </td>
   <td style="text-align:right;"> 3,127 </td>
   <td style="text-align:right;"> 426,424 </td>
   <td style="text-align:right;"> 3,617,580 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo2 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo2 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,412,814 </td>
   <td style="text-align:right;"> 80,107 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 30,124 </td>
   <td style="text-align:right;"> 211,324 </td>
   <td style="text-align:right;"> 3,148 </td>
   <td style="text-align:right;"> 249,204 </td>
   <td style="text-align:right;"> 2,986,747 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo3 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo3 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,573,998 </td>
   <td style="text-align:right;"> 84,190 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 18,210 </td>
   <td style="text-align:right;"> 249,996 </td>
   <td style="text-align:right;"> 3,579 </td>
   <td style="text-align:right;"> 320,480 </td>
   <td style="text-align:right;"> 4,250,488 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo4 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo4 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,144,340 </td>
   <td style="text-align:right;"> 71,364 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 12,009 </td>
   <td style="text-align:right;"> 304,853 </td>
   <td style="text-align:right;"> 4,605 </td>
   <td style="text-align:right;"> 336,025 </td>
   <td style="text-align:right;"> 4,873,226 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo5 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo5 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 4,991,704 </td>
   <td style="text-align:right;"> 68,566 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 11,533 </td>
   <td style="text-align:right;"> 404,111 </td>
   <td style="text-align:right;"> 6,268 </td>
   <td style="text-align:right;"> 406,501 </td>
   <td style="text-align:right;"> 5,888,730 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo6 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo6 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 3,372,007 </td>
   <td style="text-align:right;"> 41,039 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 6,472 </td>
   <td style="text-align:right;"> 323,439 </td>
   <td style="text-align:right;"> 4,908 </td>
   <td style="text-align:right;"> 270,457 </td>
   <td style="text-align:right;"> 4,018,351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo7 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo7 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,548,624 </td>
   <td style="text-align:right;"> 27,555 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 4,872 </td>
   <td style="text-align:right;"> 269,750 </td>
   <td style="text-align:right;"> 5,404 </td>
   <td style="text-align:right;"> 212,448 </td>
   <td style="text-align:right;"> 3,068,672 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polysome_786O_noVHL_noEIF4E2_g2_1_NA_ribo8 </td>
   <td style="text-align:left;"> HP5 </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> noEIF4E2 </td>
   <td style="text-align:left;"> g2 </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ribo8 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 2,377,121 </td>
   <td style="text-align:right;"> 27,535 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 6,652 </td>
   <td style="text-align:right;"> 287,538 </td>
   <td style="text-align:right;"> 3,245 </td>
   <td style="text-align:right;"> 204,979 </td>
   <td style="text-align:right;"> 2,907,083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 12,789,148 </td>
   <td style="text-align:right;"> 266,417 </td>
   <td style="text-align:right;"> 455 </td>
   <td style="text-align:right;"> 75,248 </td>
   <td style="text-align:right;"> 1,317,755 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,518,132 </td>
   <td style="text-align:right;"> 15,967,155 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 12,760,099 </td>
   <td style="text-align:right;"> 270,681 </td>
   <td style="text-align:right;"> 491 </td>
   <td style="text-align:right;"> 83,051 </td>
   <td style="text-align:right;"> 1,294,260 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,572,128 </td>
   <td style="text-align:right;"> 15,980,710 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_N_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 12,577,923 </td>
   <td style="text-align:right;"> 272,867 </td>
   <td style="text-align:right;"> 420 </td>
   <td style="text-align:right;"> 75,101 </td>
   <td style="text-align:right;"> 1,324,602 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,488,537 </td>
   <td style="text-align:right;"> 15,739,450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 14,942,365 </td>
   <td style="text-align:right;"> 309,639 </td>
   <td style="text-align:right;"> 607 </td>
   <td style="text-align:right;"> 86,114 </td>
   <td style="text-align:right;"> 1,543,982 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,765,612 </td>
   <td style="text-align:right;"> 18,648,319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 11,635,701 </td>
   <td style="text-align:right;"> 241,935 </td>
   <td style="text-align:right;"> 409 </td>
   <td style="text-align:right;"> 71,188 </td>
   <td style="text-align:right;"> 1,341,123 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,516,966 </td>
   <td style="text-align:right;"> 14,807,322 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_N_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 9,722,402 </td>
   <td style="text-align:right;"> 194,694 </td>
   <td style="text-align:right;"> 391 </td>
   <td style="text-align:right;"> 61,724 </td>
   <td style="text-align:right;"> 950,438 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,154,632 </td>
   <td style="text-align:right;"> 12,084,281 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_H_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 16,871,979 </td>
   <td style="text-align:right;"> 329,313 </td>
   <td style="text-align:right;"> 547 </td>
   <td style="text-align:right;"> 102,370 </td>
   <td style="text-align:right;"> 1,649,393 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,823,258 </td>
   <td style="text-align:right;"> 20,776,860 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_H_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 16,383,521 </td>
   <td style="text-align:right;"> 337,104 </td>
   <td style="text-align:right;"> 536 </td>
   <td style="text-align:right;"> 107,490 </td>
   <td style="text-align:right;"> 1,583,002 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,941,756 </td>
   <td style="text-align:right;"> 20,353,409 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_VHL_HIF1B_H_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 12,060,821 </td>
   <td style="text-align:right;"> 244,193 </td>
   <td style="text-align:right;"> 370 </td>
   <td style="text-align:right;"> 68,446 </td>
   <td style="text-align:right;"> 1,099,290 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,274,254 </td>
   <td style="text-align:right;"> 14,747,374 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_H_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 17,123,672 </td>
   <td style="text-align:right;"> 323,456 </td>
   <td style="text-align:right;"> 487 </td>
   <td style="text-align:right;"> 91,851 </td>
   <td style="text-align:right;"> 1,465,375 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,691,988 </td>
   <td style="text-align:right;"> 20,696,829 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_H_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 11,649,975 </td>
   <td style="text-align:right;"> 224,031 </td>
   <td style="text-align:right;"> 313 </td>
   <td style="text-align:right;"> 62,574 </td>
   <td style="text-align:right;"> 1,050,845 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,193,025 </td>
   <td style="text-align:right;"> 14,180,763 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_RCC4_noVHL_HIF1B_H_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> RCC4 </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> H </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 12,644,254 </td>
   <td style="text-align:right;"> 236,978 </td>
   <td style="text-align:right;"> 303 </td>
   <td style="text-align:right;"> 79,726 </td>
   <td style="text-align:right;"> 1,171,877 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,344,309 </td>
   <td style="text-align:right;"> 15,477,447 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_HIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6,521,700 </td>
   <td style="text-align:right;"> 133,686 </td>
   <td style="text-align:right;"> 261 </td>
   <td style="text-align:right;"> 28,606 </td>
   <td style="text-align:right;"> 942,782 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 700,294 </td>
   <td style="text-align:right;"> 8,327,329 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_HIF1B_N_2 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6,748,437 </td>
   <td style="text-align:right;"> 142,097 </td>
   <td style="text-align:right;"> 269 </td>
   <td style="text-align:right;"> 27,046 </td>
   <td style="text-align:right;"> 901,817 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 875,970 </td>
   <td style="text-align:right;"> 8,695,636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_HIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6,159,525 </td>
   <td style="text-align:right;"> 130,540 </td>
   <td style="text-align:right;"> 253 </td>
   <td style="text-align:right;"> 26,464 </td>
   <td style="text-align:right;"> 1,025,021 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 802,290 </td>
   <td style="text-align:right;"> 8,144,093 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_HIF1B_N_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 7,670,024 </td>
   <td style="text-align:right;"> 166,995 </td>
   <td style="text-align:right;"> 293 </td>
   <td style="text-align:right;"> 31,183 </td>
   <td style="text-align:right;"> 1,051,928 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 918,713 </td>
   <td style="text-align:right;"> 9,839,136 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_HIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6,612,060 </td>
   <td style="text-align:right;"> 143,839 </td>
   <td style="text-align:right;"> 213 </td>
   <td style="text-align:right;"> 26,536 </td>
   <td style="text-align:right;"> 1,061,550 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 776,444 </td>
   <td style="text-align:right;"> 8,620,642 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_HIF1B_N_2 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5,979,813 </td>
   <td style="text-align:right;"> 134,502 </td>
   <td style="text-align:right;"> 222 </td>
   <td style="text-align:right;"> 24,696 </td>
   <td style="text-align:right;"> 818,841 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 716,473 </td>
   <td style="text-align:right;"> 7,674,547 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_HIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5,310,042 </td>
   <td style="text-align:right;"> 114,411 </td>
   <td style="text-align:right;"> 180 </td>
   <td style="text-align:right;"> 28,159 </td>
   <td style="text-align:right;"> 753,813 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 616,508 </td>
   <td style="text-align:right;"> 6,823,113 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_HIF1B_N_4 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> HIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5,396,175 </td>
   <td style="text-align:right;"> 109,759 </td>
   <td style="text-align:right;"> 209 </td>
   <td style="text-align:right;"> 25,675 </td>
   <td style="text-align:right;"> 726,677 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 619,083 </td>
   <td style="text-align:right;"> 6,877,578 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_noHIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6,262,089 </td>
   <td style="text-align:right;"> 153,939 </td>
   <td style="text-align:right;"> 271 </td>
   <td style="text-align:right;"> 31,516 </td>
   <td style="text-align:right;"> 1,158,956 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 707,444 </td>
   <td style="text-align:right;"> 8,314,215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_noHIF1B_N_2 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5,782,020 </td>
   <td style="text-align:right;"> 139,342 </td>
   <td style="text-align:right;"> 223 </td>
   <td style="text-align:right;"> 21,859 </td>
   <td style="text-align:right;"> 757,813 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 638,463 </td>
   <td style="text-align:right;"> 7,339,720 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_VHL_noHIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> VHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 4,966,235 </td>
   <td style="text-align:right;"> 109,226 </td>
   <td style="text-align:right;"> 234 </td>
   <td style="text-align:right;"> 18,390 </td>
   <td style="text-align:right;"> 1,026,448 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 619,531 </td>
   <td style="text-align:right;"> 6,740,064 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_noHIF1B_N_1 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 9,449,089 </td>
   <td style="text-align:right;"> 219,167 </td>
   <td style="text-align:right;"> 361 </td>
   <td style="text-align:right;"> 46,334 </td>
   <td style="text-align:right;"> 1,718,211 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,139,727 </td>
   <td style="text-align:right;"> 12,572,889 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_noHIF1B_N_2 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 7,630,295 </td>
   <td style="text-align:right;"> 177,704 </td>
   <td style="text-align:right;"> 331 </td>
   <td style="text-align:right;"> 37,628 </td>
   <td style="text-align:right;"> 1,231,971 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 999,901 </td>
   <td style="text-align:right;"> 10,077,830 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> total_786O_noVHL_noHIF1B_N_3 </td>
   <td style="text-align:left;"> 5' end-Seq </td>
   <td style="text-align:left;"> 786O </td>
   <td style="text-align:left;"> noVHL </td>
   <td style="text-align:left;"> noHIF1B </td>
   <td style="text-align:left;"> EIF4E2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> N </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 8,779,285 </td>
   <td style="text-align:right;"> 180,624 </td>
   <td style="text-align:right;"> 347 </td>
   <td style="text-align:right;"> 37,222 </td>
   <td style="text-align:right;"> 1,410,827 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1,035,127 </td>
   <td style="text-align:right;"> 11,443,432 </td>
  </tr>
</tbody>
</table>



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
##  [9] GenomicFeatures_1.40.0      AnnotationDbi_1.50.0       
## [11] GenomicAlignments_1.24.0    Rsamtools_2.4.0            
## [13] Biostrings_2.56.0           XVector_0.28.0             
## [15] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
## [17] matrixStats_0.56.0          Biobase_2.48.0             
## [19] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
## [21] IRanges_2.22.1              S4Vectors_0.26.0           
## [23] BiocGenerics_0.34.0         rmarkdown_2.2              
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.4.1             bit64_0.9-7            viridisLite_0.3.0     
##  [4] assertthat_0.2.1       askpass_1.1            highr_0.8             
##  [7] BiocFileCache_1.12.0   blob_1.2.1             GenomeInfoDbData_1.2.3
## [10] yaml_2.2.1             progress_1.2.2         pillar_1.4.4          
## [13] RSQLite_2.2.0          lattice_0.20-41        glue_1.4.1            
## [16] digest_0.6.25          rvest_0.3.5            colorspace_1.4-1      
## [19] htmltools_0.4.0        Matrix_1.2-18          XML_3.99-0.3          
## [22] pkgconfig_2.0.3        biomaRt_2.44.0         zlibbioc_1.34.0       
## [25] purrr_0.3.4            scales_1.1.1           webshot_0.5.2         
## [28] BiocParallel_1.22.0    tibble_3.0.1           openssl_1.4.1         
## [31] farver_2.0.3           generics_0.0.2         ellipsis_0.3.1        
## [34] withr_2.2.0            crayon_1.3.4           memoise_1.1.0         
## [37] evaluate_0.14          xml2_1.3.2             tools_4.0.0           
## [40] prettyunits_1.1.1      hms_0.5.3              lifecycle_0.2.0       
## [43] munsell_0.5.0          compiler_4.0.0         rlang_0.4.6           
## [46] grid_4.0.0             RCurl_1.98-1.2         rstudioapi_0.11       
## [49] rappdirs_0.3.1         labeling_0.3           bitops_1.0-6          
## [52] gtable_0.3.0           DBI_1.1.0              curl_4.3              
## [55] R6_2.4.1               rtracklayer_1.48.0     bit_1.1-15.2          
## [58] readr_1.3.1            stringi_1.4.6          Rcpp_1.0.4.6          
## [61] vctrs_0.3.1            dbplyr_1.4.4           tidyselect_1.1.0      
## [64] xfun_0.14
```



