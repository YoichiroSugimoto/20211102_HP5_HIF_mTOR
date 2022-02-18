---
title: "sq-0 Sequence data submission"
author: "Yoichiro Sugimoto"
date: "18 February, 2022"
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



```r
## Specify the number of CPUs to be used
processors <- 4
temp <- sapply(list.files("../functions", full.names = TRUE), source)

results.dir <- file.path("../../results")

processed.fq.dir <- file.path(results.dir, "s1-processed_fastq")
processed.fq.step3.1.dir <- file.path(processed.fq.dir, "s1-1-Step3-1")

sq.dir <- file.path(results.dir, "sq-for-publication")
sq0.dir <- file.path(sq.dir, "sq0-data-submission")

create.dirs(c(
    sq.dir,
    sq0.dir
))
```

# Preparation



```r
read.process.sample.dt <- file.path(
    "../../data/sample_data/to_process_reads_sample_file.csv"
) %>%
    fread

## Exclude poor quality data due to the library pooling (i.e. HP5 of RCC4 or RCC4 VHL clone 1 with Torin 1)
read.process.sample.dt <- read.process.sample.dt[
    !grepl(
         "polysome_RCC4_(noVHL|VHL)_EIF4E2_NA_1_Torin1",
         sample_name_with_lane
     )
]

read.process.sample.dt[, `:=`(
    R1_fastq_name = paste0(sample_name_with_lane, "_R1.fastq.gz"),
    R2_fastq_name = paste0(sample_name_with_lane, "_R2.fastq.gz")
)]


print("The number of sequence files to be uploaded:")
```

```
## [1] "The number of sequence files to be uploaded:"
```

```r
addmargins(read.process.sample.dt[, table(experiment) * 2])
```

```
## experiment
## polysome    total      Sum 
##      584      128      712
```

```r
print("The number of files in s1-1-Step3-1 directory:")
```

```
## [1] "The number of files in s1-1-Step3-1 directory:"
```

```r
list.files(processed.fq.step3.1.dir) %>% length
```

```
## [1] 712
```

```r
print("Sanity check")
```

```
## [1] "Sanity check"
```

```r
list.files(processed.fq.step3.1.dir)[!(list.files(processed.fq.step3.1.dir) %in% c(read.process.sample.dt[, R1_fastq_name], read.process.sample.dt[, R2_fastq_name]))]
```

```
## character(0)
```

```r
read.process.sample.dt <- read.process.sample.dt[order(sample_name)]

total.raw.fastq.dt <- read.process.sample.dt[experiment == "total"]
hp5.raw.fastq.dt <- read.process.sample.dt[experiment == "polysome"]

fwrite(
    total.raw.fastq.dt,
    file.path(sq0.dir, "total-fastq-data.csv")
)

fwrite(
    hp5.raw.fastq.dt,
    file.path(sq0.dir, "HP5-fastq-data.csv")
)

wideFaDt <- function(raw.fastq.dt){
    m.raw.fastq.dt <- melt(
        raw.fastq.dt,
        id.vars = "no_technical_rep_sample_name",
        measure.vars = c("R1_fastq_name", "R2_fastq_name"),
        value.name = "fastq_file_name"
    )
    
    setnames(m.raw.fastq.dt, old = "no_technical_rep_sample_name", new = "sample_name")

    m.raw.fastq.dt <- m.raw.fastq.dt[order(
        sample_name, fastq_file_name
    )]
    
    m.raw.fastq.dt[, fq_file_id := paste0("fastq_", 1:.N), by = sample_name]

    d.raw.fastq.dt <- dcast(
        m.raw.fastq.dt,
        sample_name ~ fq_file_id,
        value.var = "fastq_file_name"
    )

    d.raw.fastq.dt[is.na(d.raw.fastq.dt)] <- "none"

    sample.dt <- file.path("../../data/sample_data/processed_sample_file.csv") %>%
        fread
    
    original.sample.name.order <- sample.dt[, sample_name]

    sl.sample.dt <- merge(
        sample.dt,
        d.raw.fastq.dt,
        by = "sample_name"
    )

    sl.sample.dt[
      , sample_name := factor(sample_name, levels = original.sample.name.order)
    ]

    sl.sample.dt <- sl.sample.dt[order(sample_name)]
    
    return(sl.sample.dt)
}

total.sample.dt <- wideFaDt(total.raw.fastq.dt)
hp5.sample.dt <- wideFaDt(hp5.raw.fastq.dt)

fwrite(
    total.sample.dt,
    file.path(sq0.dir, "total-sample-data.csv")
)

fwrite(
    hp5.sample.dt,
    file.path(sq0.dir, "HP5-sample-data.csv")
)
```


# Obtaining md5 for sequence data



```r
exportMd5sum <- function(raw.fastq.dt, md5sum.file){
    fq.files <- file.path(
        processed.fq.step3.1.dir,
        paste0(
            rep(raw.fastq.dt[, sample_name_with_lane], each = 2),
            c("_R1.fastq.gz", "_R2.fastq.gz")
        )
    )

    md5sum.dt <- tools::md5sum(
                            fq.files
                        ) %>%
        stack %>%
        data.table

    md5sum.dt[, ind := basename(file.path(ind))]

    md5sum.dt <- md5sum.dt[, .(ind, values)]
    
    fwrite(
        md5sum.dt, file.path(sq0.dir, md5sum.file),
        sep = "\t", col.names = FALSE
    )
    return()
}

temp <- exportMd5sum(
    raw.fastq.dt = hp5.raw.fastq.dt,
    md5sum.file = "HP5_md5sum.txt"
)

temp <- exportMd5sum(
    raw.fastq.dt = total.raw.fastq.dt,
    md5sum.file = "5Seq_md5sum.txt"
)
```

# Ftp upload of the files



```r
uploadToFtp <- function(raw.fastq.dt, ftp.url, processors, user.pwd = "aexpress:aexpress1", verbose = FALSE){

    fq.files <- file.path(
        processed.fq.step3.1.dir,
        paste0(
            rep(raw.fastq.dt[, sample_name_with_lane], each = 2),
            c("_R1.fastq.gz", "_R2.fastq.gz")
        )
    )

    if(!all(file.exists(fq.files))){
        print("could not find the following files")
        print(fq.files[!file.exists(fq.files)])
    } else {}

    if(verbose == FALSE){
        temp <- mclapply(
            fq.files,
            function(x){
                RCurl::ftpUpload(
                           what = x,
                           to = paste0(ftp.url, basename(x)),
                           userpwd = user.pwd
                       )
            },
            mc.cores = processors
        )
    } else {
        temp <- lapply(
            fq.files,
            function(x){
                print(x)
                RCurl::ftpUpload(
                           what = x,
                           to = paste0(ftp.url, basename(x)),
                           userpwd = user.pwd
                       )
            }
        )
    }
    
    return(temp)
}

cut.point <- round(nrow(hp5.raw.fastq.dt)/2)
print(paste0(cut.point * 2, " files uploaded by the first function"))
```

```
## [1] "292 files uploaded by the first function"
```

```r
ftp.url <- "ftp://aexpress:aexpress1@ftp-private-2.ebi.ac.uk/E-MTAB-10689-2lfdbpko2rskk/"

temp1 <- uploadToFtp(
    hp5.raw.fastq.dt[1:cut.point],
    ftp.url = ftp.url,
    processors = processors
)

if(nrow(data.table(stack(unlist(temp1)))[ind != "OK"]) != 0){
    stop("Error 1")
} else {"OK"}
```

```
## [1] "OK"
```

```r
temp3 <- uploadToFtp(
    hp5.raw.fastq.dt[(cut.point + 1):nrow(hp5.raw.fastq.dt)],
    ftp.url = ftp.url,
    processors = processors,
    verbose = FALSE
)

if(nrow(data.table(stack(unlist(temp3)))[ind != "OK"]) != 0){
    stop("Error 3")
} else {"OK"}
```

```
## [1] "OK"
```

```r
## temp <- uploadToFtp(
##     total.raw.fastq.dt,
##     ftp.url = "ftp://aexpress:aexpress1@ftp-private-2.ebi.ac.uk/knlecp6x-d24on2pjpl62/",
##     processors = processors
## )

## if(nrow(data.table(stack(unlist(temp)))[ind != "OK"]) != 0){
##     stop("Error 2")
## } else {"OK"}
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
##  [1] Rcpp_1.0.4.6     tidyselect_1.1.0 munsell_0.5.0    colorspace_1.4-1
##  [5] R6_2.4.1         rlang_0.4.10     tools_4.0.0      grid_4.0.0      
##  [9] gtable_0.3.0     xfun_0.14        withr_2.4.1      htmltools_0.4.0 
## [13] ellipsis_0.3.1   yaml_2.2.1       digest_0.6.25    tibble_3.0.1    
## [17] lifecycle_0.2.0  crayon_1.3.4     purrr_0.3.4      vctrs_0.3.1     
## [21] glue_1.4.1       evaluate_0.14    stringi_1.4.6    compiler_4.0.0  
## [25] pillar_1.4.4     generics_0.0.2   scales_1.1.1     pkgconfig_2.0.3
```
