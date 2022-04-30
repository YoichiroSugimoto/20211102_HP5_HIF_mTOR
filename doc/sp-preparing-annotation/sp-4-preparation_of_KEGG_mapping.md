---
title: "sp-4 KEGG and gene ID mapping"
author: "Yoichiro Sugimoto"
date: "27 April, 2022"
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
## Bioconductor database package
library("org.Hs.eg.db")
library("limma")
## Specify the number of CPUs to be used
processors <- 8

temp <- sapply(list.files("../functions", full.names = TRUE), source)

set.seed(0)
```



```r
annot.dir <- normalizePath(file.path("../../annotation/"))
kegg.dir <- file.path(annot.dir, "hg38_annotation/KEGG")

create.dirs(
    c(
        kegg.dir
    )
)
```

# Download and export KEGG id and gene id mapping



```r
kegg.link <- getGeneKEGGLinks(species = "hsa")

kegg.gene.id.dt <- AnnotationDbi::select(
         org.Hs.eg.db,
         key = kegg.link$GeneID,
         keytype = "ENTREZID",
         columns = c("ENSEMBL", "SYMBOL"),
         multiVals = "list"
         ) %>%
    data.table %>%
    {.[!duplicated(ENSEMBL)]}
```

```
## 'select()' returned many:many mapping between keys and columns
```

```r
setnames(kegg.gene.id.dt, old = "ENTREZID", new = "GeneID")

kegg.pathway.name <- getKEGGPathwayNames(species = "hsa")

kegg.data.dt <- merge(
    x = data.table(kegg.pathway.name),
    y = data.table(kegg.link),
    by = "PathwayID"
) %>%
    merge(
        kegg.gene.id.dt,
        by = "GeneID",
        allow.cartesian = TRUE
    ) 

kegg.data.dt[, `:=`(
    Description = gsub(" \\- Homo sapiens \\(human\\)", "", Description)
)]

setnames(
    kegg.data.dt,
    old = c("GeneID", "PathwayID", "Description", "ENSEMBL", "SYMBOL"),
    new = c("ENTREZID", "KEGG_id", "KEGG_description", "gene_id", "gene_name")
)

fwrite(
    kegg.data.dt,
    file.path(kegg.dir, "KEGG-gene_id-link.csv")
)
```

# KEGG Brite



```r
library("KEGGREST")

brite.tf <- "br:hsa03000"

all.brite.link.dt <- keggLink("brite", "hsa") %>%
    {data.table(stack(.))}

sl.brite.entries <- c(
    "Transcription factors" = "03000",
    "Transcription machinery" = "03021",
    "Messenger RNA biogenesis" = "03019",
    "Spliceosome" = "03041",
    "Ribosome" = "03011",
    "Ribosome biogenesis" = "03009",
    "Transfer RNA biogenesis" = "03016",
    "Translation factors" = "03012",
    "Chaperones and folding catalysts" = "03110",
    "Membrane trafficking" = "04131",
    "Ubiquitin system" = "04121",
    "Proteasome" = "03051",
    "DNA replication proteins" = "03032",
    "Chromosome and associated proteins" = "03036",
    "DNA repair and recombination proteins" = "03400",
    "Mitochondrial biogenesis" = "03029"
)

extractBriteGeneMapping <- function(brite.id, brite.name){

    brite.id.map.dt <-
        AnnotationDbi::select(
                           org.Hs.eg.db,
                           key = gsub(
                               "hsa:",
                               "",
                               all.brite.link.dt[
                                   values == paste0("br:hsa", brite.id), ind
                               ]),
                           keytype = "ENTREZID",
                           columns = c("ENSEMBL", "SYMBOL"),
                           multiVals = "list"
                       ) %>%
        data.table

    brite.id.map.dt[, BRITE_name := brite.name]

    return(brite.id.map.dt)
}

brite.gene.link.dt <- mapply(
    extractBriteGeneMapping,
    sl.brite.entries,
    names(sl.brite.entries),
    SIMPLIFY = FALSE
) %>% rbindlist
```

```
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
## 'select()' returned 1:many mapping between keys and columns
```

```r
setnames(
    brite.gene.link.dt,
    old = c("ENSEMBL"),
    new = c("gene_id")
)

fwrite(
    brite.gene.link.dt,
    file.path(kegg.dir, "KEGG_BRITE-gene_id-link.csv")
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] KEGGREST_1.28.0      knitr_1.28           stringr_1.4.0       
##  [4] magrittr_1.5         data.table_1.12.8    dplyr_1.0.0         
##  [7] khroma_1.3.0         ggplot2_3.3.1        limma_3.44.1        
## [10] org.Hs.eg.db_3.11.4  AnnotationDbi_1.50.0 IRanges_2.22.1      
## [13] S4Vectors_0.26.0     Biobase_2.48.0       BiocGenerics_0.34.0 
## [16] rmarkdown_2.2       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6      XVector_0.28.0    pillar_1.4.4      compiler_4.0.0   
##  [5] zlibbioc_1.34.0   tools_4.0.0       digest_0.6.25     bit_1.1-15.2     
##  [9] RSQLite_2.2.0     evaluate_0.14     memoise_1.1.0     lifecycle_0.2.0  
## [13] tibble_3.0.1      gtable_0.3.0      png_0.1-7         pkgconfig_2.0.3  
## [17] rlang_0.4.10      DBI_1.1.0         curl_4.3          yaml_2.2.1       
## [21] xfun_0.14         httr_1.4.2        withr_2.4.1       Biostrings_2.56.0
## [25] generics_0.0.2    vctrs_0.3.1       tidyselect_1.1.0  bit64_0.9-7      
## [29] grid_4.0.0        glue_1.4.1        R6_2.4.1          purrr_0.3.4      
## [33] blob_1.2.1        scales_1.1.1      htmltools_0.4.0   ellipsis_0.3.1   
## [37] colorspace_1.4-1  stringi_1.4.6     munsell_0.5.0     crayon_1.3.4
```
