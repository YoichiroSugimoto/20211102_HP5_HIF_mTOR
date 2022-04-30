s1-1 Preprocessing of fastq files
================
Yoichiro Sugimoto
27 April, 2022

  - [Strategy](#strategy)
  - [Set up](#set-up)
      - [Step 0. Prepare sample data](#step-0.-prepare-sample-data)
  - [Step 0. Copy raw fastq files from the data
    storage](#step-0.-copy-raw-fastq-files-from-the-data-storage)
  - [Step 1. Extract UMI from reads](#step-1.-extract-umi-from-reads)
  - [Step 2. Demultiplex reads](#step-2.-demultiplex-reads)
  - [Step 3-1. Trim adapter](#step-3-1.-trim-adapter)
  - [Delete unnecessary fq files](#delete-unnecessary-fq-files)
  - [Session Info](#session-info)

# Strategy

The fastq files will be processed before the alignments to the human
genome. The following steps will be performed:

  - Step 1: Extract UMIs from reads
  - Step 2: Demultiplex reads
  - Step 3: Process index
      - Step 3-1
          - Remove constant regions
          - Demultplex library
      - Step 3-2
          - Copy the first base of reads to read.name
  - Step 4: Map reads to rRNAs and ERCC RNAs
      - Collapse technical replicates
      - ERCC spike-in mapped reads will be used later

s1-1 will perform from Step 1 to Step 3-1, and the pre-processed files
will be submitted to ArrayExpress.

# Set up

``` r
processors <- 24

library("ShortRead")
library("stringdist") 
library("matrixStats")

temp <- sapply(list.files("../functions", full.names = TRUE), source)
```

``` r
sample.file <- file.path("../../data/sample_data/20210306_sample-data.csv")
lib2fqfilename.file <- file.path("../../data/sample_data/20210306_libname_to_fq_filename.csv")
fq.dir.file <- file.path("../../data/sample_data/20210306_fastq_dirs.csv")
tso.index.file <- file.path("../../data/sample_data/20210306_TSO-index_1.csv")
processed.sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

results.dir <- file.path("../../results")
processed.fq.dir <- file.path(results.dir, "s1-processed_fastq")

## processed.fq.log.dir <- file.path(processed.fq.dir, "log")
processed.fq.step0.dir <- file.path(processed.fq.dir, "s1-1-Step0")
processed.fq.step1.dir <- file.path(processed.fq.dir, "s1-1-Step1")
processed.fq.step2.dir <- file.path(processed.fq.dir, "s1-1-Step2") 
processed.fq.step3.1.dir <- file.path(processed.fq.dir, "s1-1-Step3-1")

index.dir <- file.path(processed.fq.dir, "index-info")

create.dirs(
    dirs = c(
        results.dir,
        processed.fq.dir,
        processed.fq.step0.dir,
        processed.fq.step1.dir,
        processed.fq.step2.dir,
        processed.fq.step3.1.dir,
        index.dir
    )
)
```

## Step 0. Prepare sample data

Processed sample data file will be used throughout.

``` r
sample.dt <- fread(sample.file)
sample.dt.cols <- colnames(sample.dt)

sample.dt[is.na(sample.dt)] <- "NA"

sample.dt[
  , sample_name := case_when(
        experiment == "total" ~ paste(
                          experiment, cell, VHL, HIF1B, oxygen, clone, replicate, sep = "_"
                      ),
        experiment == "polysome" ~ paste(
                          experiment, cell, VHL, EIF4E2, gRNA_id, clone, treatment,
                          fraction, replicate,
                          sep = "_"
                      )
    )
]

sample.dt <- sample.dt[
    order(
        experiment != "polysome",
        cell != "RCC4",
        EIF4E2 != "EIF4E2",
        oxygen != "N",
        HIF1B != "HIF1B",
        VHL != "VHL",
        treatment != "NA",
        gRNA_id,
        clone,
        replicate,
        fraction
    )
]

if(nrow(sample.dt[duplicated(sample_name)]) != 0){
    stop("sample_name duplicated")
} else {"OK"}
```

    ## [1] "OK"

``` r
sample.dt <- sample.dt[, c("sample_name", sample.dt.cols), with = FALSE]

## No technical replicate: preprocessed sample file
pre.processed.sample.dt <- copy(sample.dt[replicate == 1])
pre.processed.sample.dt[, `:=`(
    sample_name = gsub("_[[:digit:]]$", "", sample_name),
    replicate = NULL
)]

## Some data are not analysed due to the incompleteness
processed.sample.dt <- pre.processed.sample.dt[
    !grepl("_1_Torin1_ribo0", sample_name)
]

fwrite(processed.sample.dt, processed.sample.file)
```

``` r
read.process.sample.dt <- copy(sample.dt)

lib2fqfilename.dt <- fread(lib2fqfilename.file)

tso.index.dt <- fread(tso.index.file)

read.process.sample.dt <- merge(
    read.process.sample.dt,
    lib2fqfilename.dt,
    by = "library_ID",
    allow.cartesian = TRUE
) %>% {merge(
          .,
          tso.index.dt,
          by = "TSO"
      )}

## Add lane index to sample_name
read.process.sample.dt[, `:=`(
    lane_name = str_extract(R1_file_name, "SUG[[:digit:]][[:digit:]][[:digit:]]"),
    no_technical_rep_sample_name = gsub("_[[:digit:]]$", "", sample_name),
    sample_name_with_lane = paste0(
        sample_name, "_", 
        str_extract(R1_file_name, "SUG[[:digit:]][[:digit:]][[:digit:]]")
    )
)]

read.process.sample.dt <- read.process.sample.dt[
    order(
        experiment != "polysome",
        cell != "RCC4",
        EIF4E2 != "EIF4E2",
        oxygen != "N",
        HIF1B != "HIF1B",
        VHL != "VHL",
        treatment != "NA",
        gRNA_id,
        clone,
        replicate,
        fraction
    )
]

length(read.process.sample.dt[
    experiment == "total",
    unique(lane_name)
]) %>%
{print(paste0("The number of lanes used for 5' end-Seq: ", .))}
```

    ## [1] "The number of lanes used for 5' end-Seq: 3"

``` r
length(read.process.sample.dt[
    experiment == "polysome",
    unique(lane_name)
]) %>%
{print(paste0("The number of lanes used for HP5: ",.))}
```

    ## [1] "The number of lanes used for HP5: 2"

``` r
fwrite(
    read.process.sample.dt,
    file.path("../../data/sample_data/to_process_reads_sample_file.csv")
)
```

``` r
fq.dir.dt <- fread(fq.dir.file)

lib2fqfilename.dt <- lib2fqfilename.dt[
    library_ID %in% unique(sample.dt[, library_ID])
]

fq.dir.dt <- fq.dir.dt[R1_file_name %in% lib2fqfilename.dt[, R1_file_name]]

fq.dir.dt[, full_path := file.path(file_path, R1_file_name)]

print("The following files will be processed:")
```

    ## [1] "The following files will be processed:"

``` r
print(fq.dir.dt)
```

    ##                           R1_file_name
    ##  1: SUG424A5_S109_L003_R1_001.fastq.gz
    ##  2: SUG424A6_S110_L003_R1_001.fastq.gz
    ##  3: SUG424A7_S111_L003_R1_001.fastq.gz
    ##  4: SUG424A10_S99_L003_R1_001.fastq.gz
    ##  5:  SUG341A1_S23_L008_R1_001.fastq.gz
    ##  6: SUG341A2_S275_L008_R1_001.fastq.gz
    ##  7: SUG341A3_S210_L008_R1_001.fastq.gz
    ##  8: SUG341A6_S300_L008_R1_001.fastq.gz
    ##  9:  SUG341A7_S13_L001_R1_001.fastq.gz
    ## 10:  SUG341A8_S14_L001_R1_001.fastq.gz
    ## 11:  SUG341A9_S15_L001_R1_001.fastq.gz
    ## 12:  SUG341A10_S1_L001_R1_001.fastq.gz
    ## 13:  SUG341A11_S2_L001_R1_001.fastq.gz
    ## 14:  SUG341A12_S3_L001_R1_001.fastq.gz
    ## 15:  SUG341A13_S4_L001_R1_001.fastq.gz
    ## 16:  SUG341A15_S6_L001_R1_001.fastq.gz
    ## 17:  SUG341A16_S7_L001_R1_001.fastq.gz
    ## 18:  SUG341A17_S8_L001_R1_001.fastq.gz
    ## 19:  SUG341A18_S9_L001_R1_001.fastq.gz
    ## 20: SUG968A11_S68_L007_R1_001.fastq.gz
    ## 21: SUG968A12_S69_L007_R1_001.fastq.gz
    ## 22:  SUG968A1_S97_L007_R1_001.fastq.gz
    ## 23:  SUG968A2_S98_L007_R1_001.fastq.gz
    ## 24:  SUG968A3_S99_L007_R1_001.fastq.gz
    ## 25: SUG968A4_S100_L007_R1_001.fastq.gz
    ## 26: SUG968A6_S102_L007_R1_001.fastq.gz
    ## 27: SUG968A7_S103_L007_R1_001.fastq.gz
    ## 28:  SUG968A8_S71_L007_R1_001.fastq.gz
    ## 29:  SUG968A9_S72_L007_R1_001.fastq.gz
    ##                           R1_file_name
    ##                                                                                                                     file_path
    ##  1: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq
    ##  2: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq
    ##  3: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq
    ##  4: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq
    ##  5: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq
    ##  6: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq
    ##  7: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq
    ##  8: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq
    ##  9: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 10: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 11: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 12: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 13: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 14: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 15: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 16: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 17: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 18: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 19: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq
    ## 20: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq
    ## 21: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq
    ## 22: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 23: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 24: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 25: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 26: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 27: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq
    ## 28: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq
    ## 29: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq
    ##                                                                                                                     file_path
    ##                                                                                                                                                        full_path
    ##  1: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq/SUG424A5_S109_L003_R1_001.fastq.gz
    ##  2: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq/SUG424A6_S110_L003_R1_001.fastq.gz
    ##  3: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq/SUG424A7_S111_L003_R1_001.fastq.gz
    ##  4: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM18011/190806_K00102_0375_BH5KJ7BBXY/fastq/SUG424A10_S99_L003_R1_001.fastq.gz
    ##  5:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq/SUG341A1_S23_L008_R1_001.fastq.gz
    ##  6: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq/SUG341A2_S275_L008_R1_001.fastq.gz
    ##  7: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq/SUG341A3_S210_L008_R1_001.fastq.gz
    ##  8: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/190906_K00102_0391_BHC3NNBBXY/fastq/SUG341A6_S300_L008_R1_001.fastq.gz
    ##  9:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A7_S13_L001_R1_001.fastq.gz
    ## 10:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A8_S14_L001_R1_001.fastq.gz
    ## 11:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A9_S15_L001_R1_001.fastq.gz
    ## 12:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A10_S1_L001_R1_001.fastq.gz
    ## 13:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A11_S2_L001_R1_001.fastq.gz
    ## 14:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A12_S3_L001_R1_001.fastq.gz
    ## 15:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A13_S4_L001_R1_001.fastq.gz
    ## 16:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A15_S6_L001_R1_001.fastq.gz
    ## 17:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A16_S7_L001_R1_001.fastq.gz
    ## 18:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A17_S8_L001_R1_001.fastq.gz
    ## 19:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM19091/200622_K00102_0486_AHGW7LBBXY/fastq/SUG341A18_S9_L001_R1_001.fastq.gz
    ## 20: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq/SUG968A11_S68_L007_R1_001.fastq.gz
    ## 21: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq/SUG968A12_S69_L007_R1_001.fastq.gz
    ## 22:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A1_S97_L007_R1_001.fastq.gz
    ## 23:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A2_S98_L007_R1_001.fastq.gz
    ## 24:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A3_S99_L007_R1_001.fastq.gz
    ## 25: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A4_S100_L007_R1_001.fastq.gz
    ## 26: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A6_S102_L007_R1_001.fastq.gz
    ## 27: /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/201001_K00371_0406_AHJTT2BBXY/fastq/SUG968A7_S103_L007_R1_001.fastq.gz
    ## 28:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq/SUG968A8_S71_L007_R1_001.fastq.gz
    ## 29:  /camp/lab/ratcliffep/data/STPs/babs/inputs/yoichiro.sugimoto/robert.goldstone/PM20117/200915_K00102_0508_BHFH3TBBXY/fastq/SUG968A9_S72_L007_R1_001.fastq.gz
    ##                                                                                                                                                        full_path

``` r
print(lib2fqfilename.dt)
```

    ##     library_ID                       R1_file_name
    ##  1:     LYS17A SUG424A5_S109_L003_R1_001.fastq.gz
    ##  2:     LYS17B SUG424A6_S110_L003_R1_001.fastq.gz
    ##  3:     LYS17C SUG424A7_S111_L003_R1_001.fastq.gz
    ##  4:     LYS18C SUG424A10_S99_L003_R1_001.fastq.gz
    ##  5:     LYS17A  SUG341A1_S23_L008_R1_001.fastq.gz
    ##  6:     LYS17B SUG341A2_S275_L008_R1_001.fastq.gz
    ##  7:     LYS17C SUG341A3_S210_L008_R1_001.fastq.gz
    ##  8:     LYS18C SUG341A6_S300_L008_R1_001.fastq.gz
    ##  9:     LYS23A  SUG341A7_S13_L001_R1_001.fastq.gz
    ## 10:     LYS23B  SUG341A8_S14_L001_R1_001.fastq.gz
    ## 11:     LYS23C  SUG341A9_S15_L001_R1_001.fastq.gz
    ## 12:     LYS23D  SUG341A10_S1_L001_R1_001.fastq.gz
    ## 13:     LYS23E  SUG341A11_S2_L001_R1_001.fastq.gz
    ## 14:     LYS23F  SUG341A12_S3_L001_R1_001.fastq.gz
    ## 15:     LYS23G  SUG341A13_S4_L001_R1_001.fastq.gz
    ## 16:     LYS23I  SUG341A15_S6_L001_R1_001.fastq.gz
    ## 17:     LYS23J  SUG341A16_S7_L001_R1_001.fastq.gz
    ## 18:     LYS23K  SUG341A17_S8_L001_R1_001.fastq.gz
    ## 19:     LYS23L  SUG341A18_S9_L001_R1_001.fastq.gz
    ## 20:     LYS24N SUG968A11_S68_L007_R1_001.fastq.gz
    ## 21:     LYS24O SUG968A12_S69_L007_R1_001.fastq.gz
    ## 22:     LYS23J  SUG968A1_S97_L007_R1_001.fastq.gz
    ## 23:     LYS23K  SUG968A2_S98_L007_R1_001.fastq.gz
    ## 24:     LYS23L  SUG968A3_S99_L007_R1_001.fastq.gz
    ## 25:     LYS23H SUG968A4_S100_L007_R1_001.fastq.gz
    ## 26:     LYS23N SUG968A6_S102_L007_R1_001.fastq.gz
    ## 27:     LYS23O SUG968A7_S103_L007_R1_001.fastq.gz
    ## 28:     LYS23P  SUG968A8_S71_L007_R1_001.fastq.gz
    ## 29:     LYS23Q  SUG968A9_S72_L007_R1_001.fastq.gz
    ##     library_ID                       R1_file_name

# Step 0. Copy raw fastq files from the data storage

``` r
copyFastqs <- function(fq.file, processed.fq.step0.dir){
    
    out.file.full <- file.path(
        processed.fq.step0.dir,
        basename(fq.file)
    )
    
    cp.command <- paste(
        "cp", fq.file, out.file.full
    )
    system.cat(cp.command)

    cp.R2.command <- gsub("_R1_", "_R2_", cp.command)
    system.cat(cp.R2.command)

    return()
}

temp <- mclapply(
    fq.dir.dt[, full_path],
    copyFastqs,
    processed.fq.step0.dir = processed.fq.step0.dir,
    mc.cores = round(processors / 2)
)
```

# Step 1. Extract UMI from reads

``` r
step1.in.files <- file.path(
    processed.fq.step0.dir,
    unique(read.process.sample.dt[, R1_file_name])
)

print("Processing the following files for Step 1:")
```

    ## [1] "Processing the following files for Step 1:"

``` r
print(step1.in.files)
```

    ##  [1] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A8_S71_L007_R1_001.fastq.gz" 
    ##  [2] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A16_S7_L001_R1_001.fastq.gz" 
    ##  [3] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A1_S97_L007_R1_001.fastq.gz" 
    ##  [4] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A17_S8_L001_R1_001.fastq.gz" 
    ##  [5] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A2_S98_L007_R1_001.fastq.gz" 
    ##  [6] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A18_S9_L001_R1_001.fastq.gz" 
    ##  [7] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A3_S99_L007_R1_001.fastq.gz" 
    ##  [8] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A9_S72_L007_R1_001.fastq.gz" 
    ##  [9] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A6_S102_L007_R1_001.fastq.gz"
    ## [10] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A7_S103_L007_R1_001.fastq.gz"
    ## [11] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A7_S13_L001_R1_001.fastq.gz" 
    ## [12] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A8_S14_L001_R1_001.fastq.gz" 
    ## [13] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A9_S15_L001_R1_001.fastq.gz" 
    ## [14] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A10_S1_L001_R1_001.fastq.gz" 
    ## [15] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A11_S2_L001_R1_001.fastq.gz" 
    ## [16] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A12_S3_L001_R1_001.fastq.gz" 
    ## [17] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A13_S4_L001_R1_001.fastq.gz" 
    ## [18] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A4_S100_L007_R1_001.fastq.gz"
    ## [19] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A15_S6_L001_R1_001.fastq.gz" 
    ## [20] "../../results/s1-processed_fastq/s1-1-Step0/SUG424A7_S111_L003_R1_001.fastq.gz"
    ## [21] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A3_S210_L008_R1_001.fastq.gz"
    ## [22] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A11_S68_L007_R1_001.fastq.gz"
    ## [23] "../../results/s1-processed_fastq/s1-1-Step0/SUG424A10_S99_L003_R1_001.fastq.gz"
    ## [24] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A6_S300_L008_R1_001.fastq.gz"
    ## [25] "../../results/s1-processed_fastq/s1-1-Step0/SUG968A12_S69_L007_R1_001.fastq.gz"
    ## [26] "../../results/s1-processed_fastq/s1-1-Step0/SUG424A5_S109_L003_R1_001.fastq.gz"
    ## [27] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A1_S23_L008_R1_001.fastq.gz" 
    ## [28] "../../results/s1-processed_fastq/s1-1-Step0/SUG424A6_S110_L003_R1_001.fastq.gz"
    ## [29] "../../results/s1-processed_fastq/s1-1-Step0/SUG341A2_S275_L008_R1_001.fastq.gz"

``` r
extractUMI <- function(R1.fq.file, processed.fq.step1.dir){

    ## Start processing reads
    step1.in.files <- c(
        R1.fq.file,
        gsub("_R1_", "_R2_", R1.fq.file)
    )

    step1.out.files <- file.path(
        processed.fq.step1.dir,
        basename(step1.in.files)
    )

    step1.cmd <- paste0(
        "umi_tools extract",
        " -I ", step1.in.files[1],
        " --read2-in=", step1.in.files[2],
        " --bc-pattern=", "XXXXXXXXXNNNNNNN",
        " --stdout=", step1.out.files[1],
        " --read2-out=", step1.out.files[2]
    )

    step1.out <- system.cat(step1.cmd)
    ## cat(step4.out, sep = "\n")
    return()
}


temp <- mclapply(
    step1.in.files,
    extractUMI,
    processed.fq.step1.dir = processed.fq.step1.dir,
    mc.cores = processors
)
```

# Step 2. Demultiplex reads

``` r
temp.index.dt <- copy(read.process.sample.dt)

exportIndexInfo <- function(R1.fq.file, index.dir, temp.index.dt){
    temp.index1.dt <- temp.index.dt[R1_file_name == basename(R1.fq.file)]
    nrow.index1.dt <- nrow(temp.index1.dt)
    index1.dt <- rbind(
        data.table(paste0("> ", temp.index1.dt[, sample_name_with_lane])),
        data.table(paste0("^", temp.index1.dt[, index_1]))
    )[kronecker(1:nrow.index1.dt, c(0, nrow.index1.dt), "+")]

    index1.file <- file.path(index.dir, paste0(basename(R1.fq.file), ".csv"))
    fwrite(index1.dt, index1.file, col.names = FALSE)
    return()
}

temp <- sapply(
    step1.in.files,
    exportIndexInfo,
    index.dir = index.dir,
    temp.index.dt = temp.index.dt
)

demultiplexReads <- function(R1.fq.file, index.dir, processed.fq.step1.dir, processed.fq.step2.dir){
    ## Start processing reads
    step1.in.files <- c(
        R1.fq.file,
        gsub("_R1_", "_R2_", R1.fq.file)
    )

    step1.out.files <- file.path(
        processed.fq.step1.dir,
        basename(step1.in.files)
    )

    index.file <- file.path(index.dir, paste0(basename(R1.fq.file), ".csv"))
    step2.out.files <- paste0("{name}_", c("R1", "R2"), ".fastq.gz")

    step2.cmd <- paste(
        "cutadapt",
        "-g", paste0("file:", index.file),
        "-e", 0.2, # Maximum error rate
        ## "-j", processors, # Parallel computing is not supported with current arguments
        "--discard-untrimmed", # "Discard reads in which no adapter was found (cutadapt documentation)"
        "-o", file.path(processed.fq.step2.dir, step2.out.files[1]),
        "--paired-output", file.path(processed.fq.step2.dir, step2.out.files[2]),
        grep("R1", step1.out.files, value = TRUE), grep("R2", step1.out.files, value = TRUE) # input
    )

    step2.out <- system.cat(step2.cmd)
    cat(step2.out, sep = "\n")
    return()
}

temp <- mclapply(
    step1.in.files,
    demultiplexReads,
    index.dir = index.dir,
    processed.fq.step1.dir = processed.fq.step1.dir,
    processed.fq.step2.dir = processed.fq.step2.dir,
    mc.cores = processors
)
```

# Step 3-1. Trim adapter

``` r
trimAdapters <- function(analyzed.sample.name, processed.fq.step2.dir, processed.fq.step3.1.dir, processors){
    ## Start processing reads
    step3.in.files <- file.path(
        processed.fq.step2.dir,
        paste0(analyzed.sample.name, c("_R1.fastq.gz", "_R2.fastq.gz"))
    )
    
    step3.out.files <- file.path(
        processed.fq.step3.1.dir,
        basename(step3.in.files)
    )

    step3.cmd <- paste(
        "cutadapt",
        "-g", "^TTATAGGG",
        "-e", 0.2, # Maximum error rate
        "-j", processors, 
        "--discard-untrimmed", # "Discard reads in which no adapter was found (cutadapt documentation)"
        "-o", grep("_R1.fastq.gz", step3.out.files, value = TRUE),
        "--paired-output", grep("_R2.fastq.gz", step3.out.files, value = TRUE),
        grep("_R1.fastq.gz", step3.in.files, value = TRUE),
        grep("_R2.fastq.gz", step3.in.files, value = TRUE) # input
    )

    step3.out <- system.cat(step3.cmd)
    cat(step3.out, sep = "\n")
    return()
}

temp <- mclapply(
    read.process.sample.dt[, sample_name_with_lane],
    trimAdapters,
    processed.fq.step2.dir = processed.fq.step2.dir,
    processed.fq.step3.1.dir = processed.fq.step3.1.dir,
    processors = 2,
    mc.cores = round(processors / 2)
)
```

# Delete unnecessary fq files

``` r
unnecessary.samples <- pre.processed.sample.dt[, sample_name] %>%
    {.[!(. %in% processed.sample.dt[, sample_name])]}

files.to.delete <- list.files(
    processed.fq.step3.1.dir, full.names = TRUE
) %>%
    {.[grepl(
         paste0("(", paste(unnecessary.samples, collapse = "|"), ")"),
         .
     )]} 

print("The following fastq files will be deleted")
```

    ## [1] "The following fastq files will be deleted"

``` r
print(files.to.delete)
```

    ## [1] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_noVHL_EIF4E2_NA_1_Torin1_ribo0A_1_SUG968_R1.fastq.gz"
    ## [2] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_noVHL_EIF4E2_NA_1_Torin1_ribo0A_1_SUG968_R2.fastq.gz"
    ## [3] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_noVHL_EIF4E2_NA_1_Torin1_ribo0B_1_SUG968_R1.fastq.gz"
    ## [4] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_noVHL_EIF4E2_NA_1_Torin1_ribo0B_1_SUG968_R2.fastq.gz"
    ## [5] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_VHL_EIF4E2_NA_1_Torin1_ribo0A_1_SUG968_R1.fastq.gz"  
    ## [6] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_VHL_EIF4E2_NA_1_Torin1_ribo0A_1_SUG968_R2.fastq.gz"  
    ## [7] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_VHL_EIF4E2_NA_1_Torin1_ribo0B_1_SUG968_R1.fastq.gz"  
    ## [8] "../../results/s1-processed_fastq/s1-1-Step3-1/polysome_RCC4_VHL_EIF4E2_NA_1_Torin1_ribo0B_1_SUG968_R2.fastq.gz"

``` r
temp <- do.call(file.remove, list(files.to.delete))
```

# Session Info

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
    ##  [7] ggplot2_3.3.1               stringdist_0.9.5.5         
    ##  [9] ShortRead_1.46.0            GenomicAlignments_1.24.0   
    ## [11] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
    ## [13] matrixStats_0.56.0          Biobase_2.48.0             
    ## [15] Rsamtools_2.4.0             GenomicRanges_1.40.0       
    ## [17] GenomeInfoDb_1.24.0         Biostrings_2.56.0          
    ## [19] XVector_0.28.0              IRanges_2.22.1             
    ## [21] S4Vectors_0.26.0            BiocParallel_1.22.0        
    ## [23] BiocGenerics_0.34.0         rmarkdown_2.2              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0       xfun_0.14              purrr_0.3.4           
    ##  [4] lattice_0.20-41        colorspace_1.4-1       vctrs_0.3.1           
    ##  [7] generics_0.0.2         htmltools_0.4.0        yaml_2.2.1            
    ## [10] rlang_0.4.6            pillar_1.4.4           withr_2.2.0           
    ## [13] glue_1.4.1             RColorBrewer_1.1-2     jpeg_0.1-8.1          
    ## [16] GenomeInfoDbData_1.2.3 lifecycle_0.2.0        zlibbioc_1.34.0       
    ## [19] munsell_0.5.0          gtable_0.3.0           hwriter_1.3.2         
    ## [22] evaluate_0.14          latticeExtra_0.6-29    Rcpp_1.0.4.6          
    ## [25] scales_1.1.1           png_0.1-7              digest_0.6.25         
    ## [28] stringi_1.4.6          grid_4.0.0             tools_4.0.0           
    ## [31] bitops_1.0-6           RCurl_1.98-1.2         tibble_3.0.1          
    ## [34] crayon_1.3.4           pkgconfig_2.0.3        ellipsis_0.3.1        
    ## [37] Matrix_1.2-18          R6_2.4.1               compiler_4.0.0
