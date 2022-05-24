---
title: "sp-1 Preprocessing annotation"
author: "Yoichiro Sugimoto"
date: "21 May, 2022"
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

Gene annotation will be downloaded and parsed for the analyses throughout.

The following files will be downloaded:

- `GENCODE comprehensive annotation`
    - The annotation in both `gff3` and `gtf` formats will be downloaded
        - `rtracklayer` better handles the tag information in `gtf` format
        - `star` better handles `gtf` format (for at least gene counting)
- `RefSeq gtf`
  - For some genes, the annotation of `RefSeq` was more accurate than that of `GENCODE` while `GENCODE` appears to be more comprehensive in general. Thus, this script will also download RefSeq annotation and later both GENCODE and RefSeq annotations will be used for transcript assignment to TSS.
	- For example, TRIB3 or EVA1A were better annotate in `RefSeq`.

Based on these annotations, the various transcript information will be extracted including:

- `biotype` of transcripts
    - The type of transcripts such as `protein_coding`, `lncRNA`, and so on.
- `fusion gene flag`
    - The transcripts are annotated as `fusion gene` or not

This script also will select the `primary transcripts`, and export the definition in `gtf` format.

- This is a subset of `GENCODE comprehensive annotation`.
- Genes from unconventional chromosomes as well as any PAR (pseudoautosomal regions of chromosome Y) genes or genes with a PAR counterpart are exclude.


The extracted transcript annotation information is summarized and exported in the following R objects. 

- `primary.tx.dt` for GENCODE annotated transcripts
- `all.parimary.tx.dt` for both GENCODE and RefSeq annotated transcripts


# Setup



```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

### Capture `system` outputs
system <- function(...) {
    stopifnot(!any(names(list(...)) %in% "intern"))
    result <- base::system(..., intern = TRUE)
    cat(paste0(result, "\n"))
}
```


## Load packages



```r
## Bioconductor package
library("rtracklayer")
library("org.Hs.eg.db")
library("GenomicFeatures")
```


## Setup directories


```r
annot.dir <- file.path("../../annotation")
hg38.annot.dir <- file.path(
    annot.dir,
    "hg38_annotation"
)

createAnnotationDirectories <- function(folder) {
    if(dir.exists(folder) == FALSE) {
        dir.create(folder)
    }
    ## Create directories if they do not exist
    sub.directories <- file.path(folder, c("rawdata", "processed_data"))
    for(dir in sub.directories) {
        if(dir.exists(dir) == FALSE) {
            dir.create(dir)
        }
    }
}

dir.create(annot.dir)
createAnnotationDirectories(hg38.annot.dir)

annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
```


# Download and import GENCODE annotation

The [GENCODE](http://www.gencodegenes.org) gene annotation (Release 34, GRCh38.p13) is downloaded. 
The comprehensive gene annotation on the reference chromosomes are chosen (Content: Comprehensive gene annotation; Region: CHR).


## Setup data sources, filenames, and reference chromosome information



```r
gff3.url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gff3.gz"
gtf.url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz"
gencode2refseq.url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.metadata.RefSeq.gz"

file.prefix <- "GENCODE_release34_GRCh38.p13"

### reference chromosome names
conv.chromosomes <- c(as.character(seq(1,22)),"X", "Y", "M") %>%
    paste0("chr", .)
```


## Download GENCODE gene annotation



```r
gff3.file <- file.path(
    paste0(hg38.annot.dir,
           "/rawdata/",
           file.prefix,
           "_",
           format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
           ".gff3.gz")
)
download.file(
    url = gff3.url,
    destfile = gff3.file
)

gtf.file <- file.path(
    paste0(hg38.annot.dir,
           "/rawdata/",
           file.prefix,
           "_",
           format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
           ".gtf.gz")
)
download.file(
    url = gtf.url,
    destfile = gtf.file
)

## Since STAR requires unzipped file as input
ungzip.cmd <- paste(
    "gunzip", gtf.file
)

system(ungzip.cmd)
```

```r
gencode2refseq.file <- file.path(
    paste0(hg38.annot.dir,
           "/rawdata/",
           file.prefix,
           "_GENCODE2refseq_",
           format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
           ".txt.gz")
)
download.file(
    url = gencode2refseq.url,
    destfile = gencode2refseq.file
)
```


## Import GENCODE gene annotation



```r
GENCODE.annotation.gr <-import(gff3.file)
GENCODE.annotation.dt <- as.data.table(GENCODE.annotation.gr)
```


# Extract tag information from the annotation

The following tags are extracted.

- `PAR`: Used to identify duplicated genes across chromosome X and Y
- `basic`: Basic transcripts annotated by GENCODE
- MANE_Select: This is tagged for the isoforms to be used as canonical
- appris*: Additional information to rank the transcripts
- mRNA_start_NF / cds_start_NF (the mRNA / CDS start could not be confirmed): To be discarded for some of the downstream analysis
- mRNA_end_NF / cds_end_NF (the mRNA / CDS end could not be confirmed): To be fixed later
- bicistronic / readthrough_transcript: to be discarded for some of the downstream analysis

Furthermore, `refseq` tag is added to mark transcripts also annotated by `RefSeq`.



```r
## Genes in non-canonical chromosomes will be ignored
GENCODE.annotation.dt <- GENCODE.annotation.dt[seqnames %in% conv.chromosomes]

### Since the tags are aggregated in a single column as a list of vector
### the following function is used to extratc necessary tags
extractTags <- function(gff3.gr.tag.raw, tag.char){
    char.gff3.gr.tag.raw <- as.character(unlist(gff3.gr.tag.raw)) # convert to vector
    # splitted.char <- strsplit(x = char.gff3.gr.tag.raw, split = " ")
    returned.val <- if(sum(grepl(tag.char, char.gff3.gr.tag.raw)) == 0){
                        "N/A" ## In order to distinguish R's NA, N/A is used
                    } else {
                        char.gff3.gr.tag.raw[grep(tag.char, char.gff3.gr.tag.raw)]
                    }
    return(returned.val)
}

### This is a slow process; review it later
GENCODE.annotation.dt[, `:=`(
    basic_tag = sapply(tag, function(x){extractTags(x, tag.char = "^basic")}),
    PAR_tag = sapply(tag, function(x){extractTags(x, tag.char = "^PAR")}),
    MANE_Select_tag = sapply(tag, function(x){extractTags(x, tag.char = "^MANE_Select$")}),
    appris_tag = sapply(tag, function(x){extractTags(x, tag.char = "^appris_")}),
    mRNA_start_NF_tag = sapply(tag, function(x){extractTags(x, tag.char = "^mRNA_start_NF$")}),
    mRNA_end_NF_tag = sapply(tag, function(x){extractTags(x, tag.char = "^mRNA_end_NF$")}),
    cds_start_NF_tag = sapply(tag, function(x){extractTags(x, tag.char = "^cds_start_NF$")}),
    cds_end_NF_tag = sapply(tag, function(x){extractTags(x, tag.char = "^cds_end_NF$")}), 
    bicistronic_tag = sapply(tag, function(x){extractTags(x, tag.char = "^bicistronic$")}),
    readthrough_transcript_tag = sapply(tag, function(x){extractTags(x, tag.char = "^readthrough_transcript$")})
)]

par.gene.ids <- GENCODE.annotation.dt[PAR_tag == "PAR", gene_id] %>%
    unique(.)

## Genes duplicated across chromosome X and Y are removed
GENCODE.annotation.dt <- GENCODE.annotation.dt[!(gene_id %in% par.gene.ids)]

## Add refseq_tag for genes with refseq counterpart
gencode2refseq.dt <- fread(cmd = paste("zcat",  gencode2refseq.file))

setnames(
    gencode2refseq.dt,
    old = paste0("V", 1:3),
    new = c("transcript_id", "ref_nuc_id", "ref_protein_id")
)

GENCODE.annotation.dt[
  , refseq_tag := ifelse(
        transcript_id %in% gencode2refseq.dt[
                         (ref_protein_id != "") & grepl("NM_", ref_nuc_id), transcript_id
                     ],
        "refseq", "N/A"
    )
]
```

# Extract transcript / region length information to the annotation

The length of transcripts and those of transcript regions (i.e. 5' UTR, CDS, and 3'
UTR) are analyzed.


```r
### Consolidate the data.table
GENCODE.annotation.dt <- GENCODE.annotation.dt[, .(
    chromosome_name = seqnames,
    gc_strand = strand,
    gc_start = start,
    gc_end = end,
    width = width,
    type = type,
    transcript_type = transcript_type,
    gene_name = gene_name,
    gene_id = sapply(strsplit(gene_id, "\\."), "[", 1), # Gene_id's version removed
    transcript_id = transcript_id,
    phase = phase,
    exon_id = exon_id,
    exon_number = exon_number,
    basic_tag = basic_tag,
    refseq_tag = refseq_tag,
    MANE_Select_tag,
    appris_tag,
    mRNA_start_NF_tag,
    mRNA_end_NF_tag,
    cds_start_NF_tag,
    cds_end_NF_tag,
    bicistronic_tag,
    readthrough_transcript_tag
)]

### Width of transcripts and regions are calculated
rna.length.melted.dt <-
    GENCODE.annotation.dt[
      , sum(width),
        by = .(type, transcript_id)][!is.na(transcript_id)]

rna.length.dt <- dcast(
    rna.length.melted.dt,
    transcript_id ~ type,
    value.var = "V1",
    fill = 0
)

rna.length.dt <- rna.length.dt[, `:=`(
    transcript_id = transcript_id,
    tx_len = exon,
    utr5_len = five_prime_UTR,
    cds_len = CDS,
    utr3_len = three_prime_UTR
)][, .(transcript_id, tx_len, utr5_len, cds_len, utr3_len)]

setkey(GENCODE.annotation.dt, transcript_id)
setkey(rna.length.dt, transcript_id)

GENCODE.annotation.dt <- GENCODE.annotation.dt[rna.length.dt]
```


# Select primary transcripts

Primary transcripts are defined as follows:

- Non PAR genes or genes without PAR counter part
- Genes from conventional chromosomes



```r
primary.tx.dt <- GENCODE.annotation.dt[type == "transcript"]
### primary.tx.dt is consolidated
primary.tx.dt <- primary.tx.dt[, .(
    gene_name,
    gene_id,
    transcript_id,
    chromosome_name,
    biotype = "N/A", # Add later
    transcript_type,
    tx_len,
    utr5_len,
    cds_len,
    utr3_len,
    basic_tag,
    refseq_tag,
    MANE_Select_tag,
    appris_tag,
    mRNA_start_NF_tag,
    mRNA_end_NF_tag,
    cds_start_NF_tag,
    cds_end_NF_tag,
    bicistronic_tag,
    readthrough_transcript_tag
)]
```

# Define biotype of transcripts

GENCODE `transcript_type` information has too many entries, and may
not suitable for meta analysis. Thus, some of `transcript_type`
are grouped together (esp. non-coding transcripts), and the new groups
are defined as `biotype`.



```r
### Transcript_type is grouped into the following 5 classes
protein.coding <- c("protein_coding", "IG_V_gene", "TR_V_gene", "IG_C_gene",
                    "IG_J_gene", "TR_J_gene", "TR_C_gene", "IG_D_gene",
                    "TR_D_gene")
### Note that lncRNA groups are redefined later: see below
###  - tx_length > 200 nts: lncRNA
###  - tx_length <= 200 ntx: other (ncRNA)
lncRNA <- c("lncRNA", "lincRNA", "processed_transcript", "antisense","sense_intronic",
            "sense_overlapping", "3prime_overlapping_ncrna", "pseudogene",
            "transcribed_processed_pseudogene",
            "transcribed_unprocessed_pseudogene",
            "unprocessed_pseudogene", "processed_pseudogene",
            "translated_processed_pseudogene",
            "translated_unprocessed_pseudogene",
            "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA",
            "non_coding", "macro_lncRNA", "IG_pseudogene")
miRNA <- "miRNA"
rRNA <- c("rRNA", "rRNA_pseudogene")
other <- c("misc_RNA", "snRNA", "snoRNA", "retained_intron",
           "nonsense_mediated_decay", "non_stop_decay", "unitary_pseudogene",
           "polymorphic_pseudogene", "IG_V_pseudogene", "TR_V_pseudogene",
           "TR_J_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "TEC",
           "scaRNA", "transcribed_unitary_pseudogene", "scRNA", "sRNA",
           "ribozyme", "vaultRNA", "Mt_tRNA", "Mt_rRNA")
### lncRNA name dictionary is created
biotype.dict <-
    rep(c("protein_coding", "lncRNA", "miRNA", "rRNA", "other"),
        c(length(protein.coding), length(lncRNA), length(miRNA), length(rRNA), length(other))
        )
names(biotype.dict) <- c(protein.coding, lncRNA, miRNA, rRNA, other)
primary.tx.dt[, biotype := biotype.dict[transcript_type]]
### lncRNA has to be longer than 200 nts
primary.tx.dt[, biotype := ifelse((biotype == "lncRNA") & (tx_len <= 200), "other", biotype)]

## Sanity check
if(nrow(primary.tx.dt[is.na(biotype)]) != 0){
    message("Unexpected transcript type existed")
} else {"OK"}
```

```
## [1] "OK"
```

# Export primary transcript annotation in `gtf` format



```r
### Only primary transcript annotation is processed; here primary = all
GENCODE.primary.annotation.dt <- GENCODE.annotation.dt
### For IGV compatibility
GENCODE.primary.annotation.dt <- GENCODE.primary.annotation.dt[type != "transcript"]
### biotype annotation is added to GENCODE.primary.annotation.dt
txid2biotype.dict <- setNames(
    object = primary.tx.dt[, biotype],
    primary.tx.dt[, transcript_id]
    )
### To confirm the validity of using dictionary
if(
    !(class(txid2biotype.dict) == "character" &
    class(GENCODE.annotation.dt[, transcript_id]) == "character")
){stop("point 1")} else {"OK"}
```

```
## [1] "OK"
```

```r
GENCODE.primary.annotation.dt[, biotype := txid2biotype.dict[transcript_id]]

### GENCODE primary annotation is consolidated
GENCODE.primary.annotation.dt <- GENCODE.primary.annotation.dt[, .(
    type,
    chromosome_name,
    gc_start,
    gc_end,
    width,
    gc_strand,
    gene_id,
    gene_name,
    transcript_id, # Included as necessary for TxDb object generation
    phase, # Included as necessary for TxDb object generation 
    biotype,
    transcript_type,
    exon_number
)]

GENCODE.primary.annotation.dt[, exon_number := as.integer(exon_number)]
GENCODE.primary.annotation.gr <-
    makeGRangesFromDataFrame(GENCODE.primary.annotation.dt,
                             seqnames.field = "chromosome_name",
                             start.field = "gc_start",
                             end.field = "gc_end",
                             strand.field = "gc_strand",
                             keep.extra.columns = TRUE)

primary.tx.gtf <- file.path(
    paste0(hg38.annot.dir,
           "/processed_data/",
           file.prefix,
           "_primary_transcript_",
           format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
           ".gtf")
)

export(GENCODE.primary.annotation.gr, con = primary.tx.gtf, format = "gtf")
```

# Export GENCODE annotated transcript information



```r
Rdata.file <- file.path(paste0(
    hg38.annot.dir,
    "/processed_data/",
    file.prefix,
    "_primary_transcript_annotation_",
    format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
    ".rdata")
    )

save(
    conv.chromosomes, primary.tx.dt,
    file = Rdata.file
)
```


The following R objects are exported:


```r
conv.chromosomes
```

```
##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
## [19] "chr19" "chr20" "chr21" "chr22" "chrX"  "chrY"  "chrM"
```

```r
primary.tx.dt
```

```
##          gene_name         gene_id      transcript_id chromosome_name
##      1:       ARF5 ENSG00000004059 ENST00000000233.10            chr7
##      2:       M6PR ENSG00000003056  ENST00000000412.8           chr12
##      3:      ESRRA ENSG00000173153 ENST00000000442.11           chr11
##      4:      FKBP4 ENSG00000004478  ENST00000001008.6           chr12
##      5:    CYP26B1 ENSG00000003137  ENST00000001146.6            chr2
##     ---                                                              
## 227722:     DNAJC7 ENSG00000168259  ENST00000674504.1           chr17
## 227723:    HSPA12A ENSG00000165868  ENST00000674505.1           chr10
## 227724: AL354833.1 ENSG00000288598  ENST00000674506.1           chr13
## 227725:      GFPT1 ENSG00000198380  ENST00000674507.1            chr2
## 227726:      GNL3L ENSG00000130119  ENST00000674508.1            chrX
##                biotype         transcript_type tx_len utr5_len cds_len utr3_len
##      1: protein_coding          protein_coding   1032       88     543      401
##      2: protein_coding          protein_coding   2450      159     834     1457
##      3: protein_coding          protein_coding   2274      225    1272      777
##      4: protein_coding          protein_coding   3715      170    1380     2165
##      5: protein_coding          protein_coding   4732      204    1539     2989
##     ---                                                                        
## 227722:          other nonsense_mediated_decay   2854        9     120     2725
## 227723:          other nonsense_mediated_decay   5676      354     549     4773
## 227724:         lncRNA                  lncRNA   1017        0       0        0
## 227725: protein_coding          protein_coding   5906      117    1878     3911
## 227726: protein_coding          protein_coding    278      261      17        0
##         basic_tag refseq_tag MANE_Select_tag         appris_tag
##      1:     basic     refseq     MANE_Select appris_principal_1
##      2:     basic     refseq     MANE_Select appris_principal_1
##      3:     basic     refseq     MANE_Select appris_principal_1
##      4:     basic     refseq     MANE_Select appris_principal_1
##      5:     basic     refseq             N/A appris_principal_1
##     ---                                                        
## 227722:       N/A        N/A             N/A                N/A
## 227723:       N/A        N/A             N/A                N/A
## 227724:       N/A        N/A             N/A                N/A
## 227725:     basic        N/A             N/A                N/A
## 227726:       N/A        N/A             N/A                N/A
##         mRNA_start_NF_tag mRNA_end_NF_tag cds_start_NF_tag cds_end_NF_tag
##      1:               N/A             N/A              N/A            N/A
##      2:               N/A             N/A              N/A            N/A
##      3:               N/A             N/A              N/A            N/A
##      4:               N/A             N/A              N/A            N/A
##      5:               N/A             N/A              N/A            N/A
##     ---                                                                  
## 227722:               N/A             N/A              N/A            N/A
## 227723:               N/A             N/A              N/A            N/A
## 227724:               N/A             N/A              N/A            N/A
## 227725:               N/A             N/A              N/A            N/A
## 227726:               N/A     mRNA_end_NF              N/A     cds_end_NF
##         bicistronic_tag readthrough_transcript_tag
##      1:             N/A                        N/A
##      2:             N/A                        N/A
##      3:             N/A                        N/A
##      4:             N/A                        N/A
##      5:             N/A                        N/A
##     ---                                           
## 227722:             N/A                        N/A
## 227723:             N/A                        N/A
## 227724:             N/A                        N/A
## 227725:             N/A                        N/A
## 227726:             N/A                        N/A
```


# Download and process RefSeq annotation



```r
refseq.file.prefix <- "RefSeq_GRCh38.p13"

refseq.gtf.file <- list.files(
    file.path(hg38.annot.dir, "rawdata"),
    pattern = paste0(refseq.file.prefix, ".*\\.gtf\\.gz$"),
    full.names = TRUE
) 

if(length(refseq.gtf.file) == 0){    
    refseq.gtf.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20200522/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz"
    refseq.gtf.file <- file.path(
        paste0(hg38.annot.dir,
               "/rawdata/",
               refseq.file.prefix,
               "_",
               format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
               ".gtf.gz")
    )
    download.file(
        url = refseq.gtf.url,
        destfile = refseq.gtf.file
    )
} else {"RefSeq gtf already downaloded"}

refseq.tx.gr <- import(refseq.gtf.file)
refseq.tx.dt <- refseq.tx.gr %>% as.data.frame %>% data.table

## Fix chromosome names in RefSeq to conventional naming
refseq.conv.chr <- unique(refseq.tx.dt[, seqnames]) %>%
    {grep("^NC_", x = ., value = TRUE)}

temp.refseq.conv.chr2standard.chr <- str_split_fixed(refseq.conv.chr, "_", n = 2)[, 2] %>%
    {str_split_fixed(string = ., pattern = "\\.", n = 2)[, 1]} %>%
    as.integer %>%
    {setNames(
        object = paste0("chr", .), 
        nm = refseq.conv.chr
     )}

refseq.conv.chr2standard.chr <- temp.refseq.conv.chr2standard.chr %>% {case_when(
    . == "chr23" ~ "chrX",
    . == "chr24" ~ "chrY",
    . == "chr12920" ~ "chrM",
    TRUE ~ .
)}
names(refseq.conv.chr2standard.chr) <- names(temp.refseq.conv.chr2standard.chr)
refseq.tx.dt[, seqnames := refseq.conv.chr2standard.chr[seqnames]]

## Discard unnecessary annotations (since this is for tx assignment to use as secondary annotation in addition to GENCODE annotation, I discarded any unusual annotations)
sl.refseq.tx.dt <- refseq.tx.dt[
    type != "gene" &
    seqnames %in% refseq.conv.chr2standard.chr &
    grepl("(^NM_|^NR_)", transcript_id)
]

## Change gene_id to ensembl gene id
sl.refseq.tx.dt <-
    mapIds(
        org.Hs.eg.db,
        keys = gsub("^GeneID:", "", sl.refseq.tx.dt[, db_xref]),
        keytype = "ENTREZID", column = "ENSEMBL",
        multiVals = "first"
    ) %>%
    {sl.refseq.tx.dt[, gene_id := .]} %>%
    {.[!is.na(gene_id)]}
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
## Add gene name
sl.refseq.tx.dt <- inner_join(
    sl.refseq.tx.dt,
    primary.tx.dt[!duplicated(gene_id), .(gene_id, gene_name)],
    by = "gene_id"
)

## Remove columns with only NA
## https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
sl.refseq.tx.dt <- sl.refseq.tx.dt[
   ,
    which(unlist(lapply(sl.refseq.tx.dt, function(x)!all(is.na(x))))),
    with = FALSE
]

## Prapre primary tx.dt for refseq data
## Calculate tx length
refseq.primary.txdb <- makeTxDbFromGRanges(
    makeGRangesFromDataFrame(sl.refseq.tx.dt, keep.extra.columns = TRUE)
)
```

```
## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): The "phase" metadata column contains non-NA values for features of type
##   stop_codon. This information was ignored.
```

```
## Warning in .find_exon_cds(exons, cds): The following transcripts have exons that contain more than one CDS
##   (only the first CDS was kept for each exon): NM_001134939.1,
##   NM_001172437.2, NM_001184961.1, NM_001301020.1, NM_001301302.1,
##   NM_001301371.1, NM_002537.3, NM_004152.3, NM_015068.3, NM_016178.2
```

```
## Warning in .reject_transcripts(bad_tx, because): The following transcripts were rejected because they have incompatible
##   CDS and stop codons: NM_001172437.2, NM_001184961.1, NM_015068.3
```

```r
refseq.primary.tx.len.dt <- transcriptLengths(
    refseq.primary.txdb,
    with.cds_len = TRUE,
    with.utr5_len = TRUE,
    with.utr3_len = TRUE
) %>% data.table %>%
{.[, `:=`(
     tx_id = NULL,
     nexon = NULL
 )]}

setnames(refseq.primary.tx.len.dt, old = "tx_name", "transcript_id")

refseq.primary.tx.dt <- sl.refseq.tx.dt[
  , .(gene_id, gene_name, transcript_id, seqnames, gbkey)] %>%
    {.[!duplicated(.)][gbkey != "CDS"]} %>%
    merge(        
        refseq.primary.tx.len.dt,
        by = c("transcript_id", "gene_id")
    ) %>%
    {.[, `:=`(
         transcript_type = case_when(
             gbkey == "mRNA" ~ "protein_coding",
             gbkey == "ncRNA" & tx_len > 200 ~ "lncRNA",
             TRUE ~ gbkey
         ) %>%
             factor(
                 levels = c(
                     "protein_coding", "lncRNA", "ncRNA",
                     "misc_RNA", "precursor_RNA", "rRNA"
                 ), ordered = TRUE),
         basic_tag = "N/A",
         refseq_tag = "original_refseq",
         MANE_Select_tag = "N/A",
         appris_tag = "N/A",
         mRNA_start_NF_tag = "N/A",
         mRNA_end_NF_tag = "N/A",
         cds_start_NF_tag = "N/A",
         cds_end_NF_tag = "N/A",
         bicistronic_tag = "N/A",
         readthrough_transcript_tag = "N/A"
     )]} %>%
    {.[, gbkey := NULL]}


refseq.biotype.dt <- refseq.primary.tx.dt[
    order(transcript_type)][, head(.SD, 1), by = gene_id] %>%
    {.[, biotype := transcript_type]}

refseq.primary.tx.dt <- merge(
    refseq.primary.tx.dt,
    refseq.biotype.dt[!duplicated(gene_id), .(gene_id, biotype)],
    by = "gene_id"
)

setnames(refseq.primary.tx.dt, old = "seqnames", new = "chromosome_name")

all.primary.tx.dt <- rbind(
    primary.tx.dt,
    refseq.primary.tx.dt
)

## Fusion gene can confound the analysis and here I flag them
all.primary.tx.dt[, fusion_gene_flag := ifelse(
                        transcript_type == "protein_coding" &
                        grepl("-[[:alpha:]]", gene_name) &
                        !grepl("-AS[[:digit:]]$", gene_name) &
                        !grepl("^MT-", gene_name) &
                        !grepl("^HLA-", gene_name) |
                        readthrough_transcript_tag == "readthrough_transcript" |
                        bicistronic_tag == "bicistronic",
                        TRUE, FALSE
                    )]


## Add transcript_type and biotype information for refseq.gtf
sl.refseq.tx.dt <- inner_join(
    sl.refseq.tx.dt,
    all.primary.tx.dt[, .(transcript_id, transcript_type, biotype)],
    by = "transcript_id"
)

## Remove columns which are not necessary
sl.refseq.tx.dt <- sl.refseq.tx.dt[, c(
    "seqnames", "start", "end", "width", "strand", "source", "type",
    "phase", "gene_id", "gene_name", "transcript_id", "biotype", "transcript_type", "exon_number"
)]

processed.refseq.gtf <- file.path(
    annot.ps.dir,
    paste0(basename(refseq.gtf.file)) %>%
    {gsub("\\.gtf\\.gz$", "_processed.gtf", .)}
)

rtracklayer::export(
                 makeGRangesFromDataFrame(sl.refseq.tx.dt, keep.extra.columns = TRUE),
                 con = processed.refseq.gtf
             )
```


# Export GENCODE and RefSeq annotated transcript information


```r
all.primary.tx.info.file <- file.path(
    annot.ps.dir,
    "all_GENCODE_RefSeq_transcript_info.csv"
)

fwrite(all.primary.tx.dt, all.primary.tx.info.file)
```


# Create gtf file for all transcript



```r
primary.tx.gr <- rtracklayer::import(primary.tx.gtf)
refseq.tx.gr <- rtracklayer::import(processed.refseq.gtf)

all.tx.gr <- c(primary.tx.gr, refseq.tx.gr)

all.primary.tx.dt[
  , pc_gene_flag := gene_id %in% all.primary.tx.dt[biotype == "protein_coding", gene_id]
]

sl.tx.ids <- all.primary.tx.dt[
    refseq_tag != "refseq" &
    ## pc_gene_flag == (transcript_type == "protein_coding") &
    mRNA_start_NF_tag == "N/A" &
    ## mRNA_end_NF_tag == "N/A" &
    ## cds_start_NF_tag == "N/A" &
    ## cds_end_NF_tag == "N/A" &
    bicistronic_tag == "N/A" &
    readthrough_transcript_tag == "N/A" &
    fusion_gene_flag == FALSE,
    transcript_id
]

all.tx.gr <- all.tx.gr[
    mcols(all.tx.gr)[["transcript_id"]] %in% sl.tx.ids
]

all.tx.gtf <- file.path(
    annot.ps.dir,
    "all-transcript.gtf"
)

rtracklayer::export(
                 all.tx.gr,
                 con = all.tx.gtf
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.40.0 org.Hs.eg.db_3.11.4    AnnotationDbi_1.50.0  
##  [4] Biobase_2.48.0         rtracklayer_1.48.0     GenomicRanges_1.40.0  
##  [7] GenomeInfoDb_1.24.0    IRanges_2.22.1         S4Vectors_0.26.0      
## [10] BiocGenerics_0.34.0    knitr_1.28             stringr_1.4.0         
## [13] magrittr_1.5           data.table_1.12.8      dplyr_1.0.0           
## [16] khroma_1.3.0           ggplot2_3.3.1          rmarkdown_2.2         
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.4.2                  bit64_0.9-7                
##  [3] assertthat_0.2.1            askpass_1.1                
##  [5] BiocFileCache_1.12.0        blob_1.2.1                 
##  [7] GenomeInfoDbData_1.2.3      Rsamtools_2.4.0            
##  [9] yaml_2.2.1                  progress_1.2.2             
## [11] pillar_1.4.4                RSQLite_2.2.0              
## [13] lattice_0.20-41             glue_1.4.1                 
## [15] digest_0.6.25               XVector_0.28.0             
## [17] colorspace_1.4-1            htmltools_0.4.0            
## [19] Matrix_1.2-18               XML_3.99-0.3               
## [21] pkgconfig_2.0.3             biomaRt_2.44.0             
## [23] zlibbioc_1.34.0             purrr_0.3.4                
## [25] scales_1.1.1                BiocParallel_1.22.0        
## [27] tibble_3.0.1                openssl_1.4.1              
## [29] generics_0.0.2              ellipsis_0.3.1             
## [31] withr_2.4.1                 SummarizedExperiment_1.18.1
## [33] crayon_1.3.4                memoise_1.1.0              
## [35] evaluate_0.14               tools_4.0.0                
## [37] prettyunits_1.1.1           hms_0.5.3                  
## [39] lifecycle_0.2.0             matrixStats_0.56.0         
## [41] munsell_0.5.0               DelayedArray_0.14.0        
## [43] Biostrings_2.56.0           compiler_4.0.0             
## [45] rlang_0.4.10                grid_4.0.0                 
## [47] RCurl_1.98-1.2              rappdirs_0.3.1             
## [49] bitops_1.0-6                gtable_0.3.0               
## [51] curl_4.3                    DBI_1.1.0                  
## [53] R6_2.4.1                    GenomicAlignments_1.24.0   
## [55] bit_1.1-15.2                stringi_1.4.6              
## [57] Rcpp_1.0.4.6                vctrs_0.3.1                
## [59] dbplyr_1.4.4                tidyselect_1.1.0           
## [61] xfun_0.14
```
