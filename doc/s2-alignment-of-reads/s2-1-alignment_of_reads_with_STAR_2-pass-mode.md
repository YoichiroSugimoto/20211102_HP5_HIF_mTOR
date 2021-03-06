s2-1 Alignment of reads with STAR (2 pass mode)
================
Yoichiro Sugimoto
22 May, 2022

  - [Strategy](#strategy)
  - [Set up](#set-up)
  - [Align reads with STAR](#align-reads-with-star)
  - [Session Information](#session-information)

# Strategy

Sequence reads will be aligned to the human genome using STAR

# Set up

``` r
processors <- 24

temp <- sapply(list.files("../functions", full.names = TRUE), source)
```

``` r
sample.file <- file.path("../../data/sample_data/processed_sample_file.csv")

results.dir <- file.path("../../results")
processed.fq.dir <- file.path(results.dir, "s1-processed_fastq")

## Input fq files
processed.fq.step4.dir <- file.path(processed.fq.dir, "s1-1-Step4") 

annot.dir <- file.path("../../annotation")
annot.ps.dir <- file.path(annot.dir, "hg38_annotation/processed_data/")
all.tx.gtf <- file.path(
    annot.ps.dir,
    "all-transcript.gtf"
)

star.index.dir <- file.path(annot.dir, "hg38_annotation/star_indices/")
s2.alignment.dir <- file.path(results.dir, "s2-read-alignment")
star.aligned.read.dir <- file.path(s2.alignment.dir, "s2-1-a-star-aligned_reads")
star.aligned.bam.dir <- file.path(s2.alignment.dir, "s2-1-b-star-aligned_bam")

create.dirs(c(
    s2.alignment.dir,
    star.aligned.read.dir,
    star.aligned.bam.dir
))
```

# Align reads with STAR

``` r
starCommand <- function(fastq.files, start.index.dir, out.dir, sample.name, processors, additional.args = ""){

    cmd <- paste("STAR",
                 "--outFilterType", "BySJout", ## ENCODE option; reduce the number of spurious junctions
                 ## Below are ENCODE options
                 ## "--alignSJoverhangMin", 8,
                 ## "--alignSJDBoverhangMin", 1, 
                 ## "--outFilterMismatchNmax", 999,
                 ## "--outFilterMismatchNoverReadLmax", 0.04,
                 ## "--alignIntronMin", 20,
                 ## "--alignIntronMax", 1000000,
                 ## "--alignMatesGapMax", 1000000,
                 ## ENCODE option untill here
                 "--outFilterMultimapNmax", 1, # Report only uniquely mapped reads
                 ## "--chimSegmentMin", 2,
                 ## "--outFilterMismatchNmax", 3,
                 ## "--alignEndsType", "EndToEnd", ## No softclip
                 "--runThreadN", processors,
                 "--readFilesCommand", "zcat",
                 "--limitBAMsortRAM 60000000000", # from https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/small-rna/small-rna-align/resources/usr/bin/srna_align.sh
                 "--genomeDir", star.index.dir,
                 "--readFilesIn", paste(fastq.files, collapse = " "),
                 "--quantMode", "GeneCounts",
                 "--outFileNamePrefix", paste0(
                                            out.dir, "/",
                                            sample.name, "_"
                                        ),
                 additional.args
                 )
    
    return(cmd)
}



alignReadsWithSTAR <- function(sample.name,
                               processed.fq.step4.dir,
                               star.aligned.read.dir,
                               star.aligned.bam.dir,
                               star.index.dir,
                               all.tx.gtf,
                               second.run.flag,
                               counter,
                               inserted.genome.dir,
                               processors
                               ){

    ## Input
    fq.files <- file.path(
        processed.fq.step4.dir,
        paste0(sample.name, ".fastq.", 1:2, ".gz")
    )
    ## Output
    out.dir <- file.path(star.aligned.read.dir, sample.name)
    create.dir(out.dir)

    ## Bam file name
    out.bam <- file.path(
        star.aligned.bam.dir,
        paste0(sample.name, ".bam")
    )
    out.bai <- gsub(".bam$", ".bai", out.bam)

    ## If this is the first run of the second round, insert junction into the genome
    if((second.run.flag == TRUE) & (counter == 1)){
        additional.args <- paste(
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--sjdbFileChrStartEnd", paste(
                                         list.files(
                                             star.aligned.read.dir,
                                             recursive = TRUE,
                                             pattern = "SJ.out.tab$",
                                             full.names = TRUE)
                                     ),
            "--sjdbGTFfile", all.tx.gtf,
            "--sjdbOverhang", 99, # max(ReadLength)-1
            "--sjdbInsertSave", "All"
        )
        genome.dir <- star.index.dir 
    } else if ((second.run.flag == TRUE) & (counter > 1)) {
        additional.args <- paste(
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--genomeLoad LoadAndKeep"
        )
        genome.dir <-  inserted.genome.dir
    } else if (second.run.flag == FALSE){
        additional.args <- paste(
            "--outSAMtype", "None",
            "--genomeLoad LoadAndKeep"
        )
        genome.dir <- star.index.dir
    } else {
        stop("This should not happen...")
    }
    
    star.cmd <- starCommand(
        fastq.files = c(fq.files),
        start.index.dir = start.index.dir,
        out.dir = out.dir,
        sample.name = sample.name,
        processor = processors,
        additional.args = additional.args
    )
    
    system.cat(star.cmd)

    ## Move genomic coordinate mapped data to output folder
    move.cmd <- paste(
        "mv",
        file.path(out.dir, paste0(sample.name, "_Aligned.sortedByCoord.out.bam")),
        out.bam
    )

    bai.cmd <- paste(
        "samtools index",
        "-@", processors,
        out.bam,
        out.bai
    )

    if(second.run.flag == TRUE){
        system.cat(move.cmd)
        system.cat(bai.cmd)
    } else {
        
    }

    inserted.genome.dir <- file.path(out.dir, paste0(sample.name, "__STARgenome"))
    
    return(inserted.genome.dir)
}


sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]

## Clean genome from shared memory
temp.dir <- file.path(star.aligned.read.dir, "temp/")
create.dir(temp.dir)
clear.genome.cmd <- paste(
    "STAR",
    "--genomeLoad Remove",
    "--genomeDir", star.index.dir,
    "--runThreadN", processors,
    "--outFileNamePrefix", temp.dir
)

system.cat(clear.genome.cmd)
```

    ## Warning in system(cmd, intern = TRUE): running command 'STAR --genomeLoad Remove
    ## --genomeDir ../../annotation/hg38_annotation/star_indices/ --runThreadN 24 --
    ## outFileNamePrefix ../../results/s2-read-alignment/s2-1-a-star-aligned_reads/
    ## temp/' had status 105

    ## [1] "May 22 06:54:33 ..... started STAR run"
    ## [2] "May 22 06:54:33 ..... loading genome"  
    ## attr(,"status")
    ## [1] 105

``` r
## 1st pass

for(i in 1:length(sample.names)){

    sample.name <- sample.names[i]
    
    temp.inserted.genome.dir <- alignReadsWithSTAR(
        sample.name = sample.name,
        processed.fq.step4.dir = processed.fq.step4.dir,
        star.aligned.read.dir = star.aligned.read.dir,
        star.aligned.bam.dir = star.aligned.bam.dir,
        star.index.dir = star.index.dir,
        all.tx.gtf = all.tx.gtf,
        second.run.flag = FALSE,
        counter = i,
        inserted.genome.dir = "",
        processors = processors
    )

}

## 2nd pass
for(i in 1:length(sample.names)){

    sample.name <- sample.names[i]

    if(i == 1){
        inserted.genome.dir <- ""
        system.cat(clear.genome.cmd)
    } else {}

    temp.inserted.genome.dir <- alignReadsWithSTAR(
        sample.name = sample.name,
        processed.fq.step4.dir = processed.fq.step4.dir,
        star.aligned.read.dir = star.aligned.read.dir,
        star.aligned.bam.dir = star.aligned.bam.dir,
        star.index.dir = star.index.dir,
        all.tx.gtf = all.tx.gtf,
        second.run.flag = TRUE,
        counter = i,
        inserted.genome.dir = inserted.genome.dir,
        processors = processors
    )
    
    if(i == 1){

        inserted.genome.dir <- temp.inserted.genome.dir
        system.cat(clear.genome.cmd)
        
    } else {
        ## 
    }
}
```

    ## Warning in system(cmd, intern = TRUE): running command 'STAR --genomeLoad Remove
    ## --genomeDir ../../annotation/hg38_annotation/star_indices/ --runThreadN 24 --
    ## outFileNamePrefix ../../results/s2-read-alignment/s2-1-a-star-aligned_reads/
    ## temp/' had status 105

``` r
## Delete generated genome index
## unlink(inserted.genome.dir, recursive = TRUE)

## Remove genome from shared memory
system.cat(clear.genome.cmd) 
```

    ## [1] "May 22 10:41:56 ..... started STAR run"
    ## [2] "May 22 10:41:56 ..... loading genome"

``` r
unlink(temp.dir, recursive = TRUE)
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.28        stringr_1.4.0     magrittr_1.5      data.table_1.12.8
    ## [5] dplyr_1.0.0       khroma_1.3.0      ggplot2_3.3.1     rmarkdown_2.2    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6     tidyselect_1.1.0 munsell_0.5.0    colorspace_1.4-1
    ##  [5] R6_2.4.1         rlang_0.4.6      tools_4.0.0      grid_4.0.0      
    ##  [9] gtable_0.3.0     xfun_0.14        withr_2.2.0      htmltools_0.4.0 
    ## [13] ellipsis_0.3.1   yaml_2.2.1       digest_0.6.25    tibble_3.0.1    
    ## [17] lifecycle_0.2.0  crayon_1.3.4     purrr_0.3.4      vctrs_0.3.1     
    ## [21] glue_1.4.1       evaluate_0.14    stringi_1.4.6    compiler_4.0.0  
    ## [25] pillar_1.4.4     generics_0.0.2   scales_1.1.1     pkgconfig_2.0.3
