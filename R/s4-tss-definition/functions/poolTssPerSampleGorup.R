suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(
        c("-s", "--sample_group"),
        action="store",
        type = "character",        
        help = "sample group name"),    
    make_option(
        c("-b", "--bam_file_dir"),
        action="store",
        type = "character",        
        help = "bam file path"),
    make_option(
        c("--with_g_output_directory"),
        action="store",
        type = "character",
        help = "directory for TSS with untemplated G output files"),
    make_option(
        c("--without_g_output_directory"),
        action="store",
        type = "character",
        help = "directory for TSS without untemplated G output files"),
    make_option(
        c("--tss_table_output_directory"),
        action="store",
        type = "character",
        help = "directory for TSS table file"),    
    make_option(
        c("-p", "--processors"),
        action="store",
        type = "character",        
        help = "the number of CPUs to be used")    
)

opt <- parse_args(OptionParser(option_list=option_list))

poolTssPerSampleGorup <- function(sample.group, tss.bam.dir, s4.1.5.1.tss.bam.with.untemplated.G.dir, s4.1.5.2.tss.bam.without.untemplated.G.dir, s4.1.5.3.tss.table.dir, processors){

    ## Data analysis
    suppressPackageStartupMessages(require("data.table"))
    suppressPackageStartupMessages(require("magrittr"))
    suppressPackageStartupMessages(require("stringr"))
    ## Bioconductor
    suppressPackageStartupMessages(require("rtracklayer"))
    suppressPackageStartupMessages(require("GenomicAlignments"))
    ## Parallelization
    suppressPackageStartupMessages(require("parallel"))
    
    bam.files <- list.files(
        tss.bam.dir,  
        pattern = paste0(sample.group, ".*\\.bam$"), 
        full.names = TRUE  
    )
    
    readTssBam <- function(bam.file){
        param <- ScanBamParam(
            flag = scanBamFlag(isSecondaryAlignment = FALSE),
            what = c("flag"),
            tag = c("XG") # this tag will define reads containing untemplated G
        ) # uniquely mapped

        bam <- readGAlignments(bam.file, use.names = TRUE, param = param)

        utg.bam.con <- file.path(
            s4.1.5.1.tss.bam.with.untemplated.G.dir,
            basename(bam.file)
        )
        nog.bam.con <- file.path(
            s4.1.5.2.tss.bam.without.untemplated.G.dir,
            basename(bam.file)
        )

        export(bam[mcols(bam)$XG == 1], con = utg.bam.con, format = "bam")
        export(bam[mcols(bam)$XG == 0], con = nog.bam.con, format = "bam")
        
        bam.dt <- as.data.table(bam)
        bam.dt <- bam.dt[, .(seqnames, strand, start)]
        setnames(bam.dt, old = "seqnames", new = "chr")
        setkey(bam.dt, chr, strand, start)
        return(bam.dt)
    }

    tsss <- mclapply(
        bam.files,
        readTssBam,
        mc.cores = processors
    )

    tss.dt <- rbindlist(tsss)
    setkey(tss.dt, chr, strand, start)

    ## tss.count.dt <- tss.dt[, .N, by = list(chr, strand, start)]
    ## setnames(tss.count.dt, old = "N", new = "count")

    fwrite(
        tss.dt,
        file.path(s4.1.5.3.tss.table.dir, paste0(sample.group, ".csv"))
    )

    return(tss.dt)
}

temp <- poolTssPerSampleGorup(
    sample.group = opt$sample_group,
    tss.bam.dir = opt$bam_file_dir,
    s4.1.5.1.tss.bam.with.untemplated.G.dir = opt$with_g_output_directory,
    s4.1.5.2.tss.bam.without.untemplated.G.dir = opt$without_g_output_directory,
    s4.1.5.3.tss.table.dir = opt$tss_table_output_directory,
    processors = opt$processors
)
