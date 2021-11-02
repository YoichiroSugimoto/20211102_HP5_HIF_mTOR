## This function trim soft clipped bases because of the reverse transcription of mRNA cap from bam files

suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
    make_option(
        c("-s", "--sample_name"),
        action="store",
        type = "character",        
        help = "sample name"),
    make_option(
        c("-b", "--star_bam_dir"),
        action="store",
        type = "character",        
        help = "sample name"),
    make_option(
        c("-d", "--dedup_tss_bam"),
        action="store",
        type = "character",        
        help = "sample name"),    
    make_option(
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("rtracklayer"))


dedupBam <- function(sample.name, star.aligned.bam.dir, s2.2.2.dedup.tss.bam.dir, s2.2.3.dedup.bam.dir){

    input.bam <- list.files(
        star.aligned.bam.dir,
        pattern = glob2rx(paste0(sample.name, "*.bam")),
        full.names = TRUE
    )
        
    ref.bam <- list.files(
        s2.2.2.dedup.tss.bam.dir,
        pattern = glob2rx(paste0(sample.name, "*.bam")),
        full.names = TRUE
    )

    output.bam <- file.path(
        s2.2.3.dedup.bam.dir,
        gsub(".bam", ".dedup.bam", basename(input.bam))
    )

    dedup.qnames <- names(readGAlignments(ref.bam, use.names = TRUE))

    param <- ScanBamParam(
        what = c("flag", "isize", "seq", "qual")
    )

    bam <- readGAlignmentPairs(input.bam, use.names = TRUE, param=param)

    dedup.bam <- bam[names(bam) %in% dedup.qnames]

    rm(bam)
    gc()
    
    export(dedup.bam, con = output.bam, format = "bam")

    rm(dedup.bam)
    gc()
    
    return()
}

temp <- dedupBam(
    sample.name = opt$sample_name,
    star.aligned.bam.dir = opt$star_bam_dir,
    s2.2.2.dedup.tss.bam.dir = opt$dedup_tss_bam,
    s2.2.3.dedup.bam.dir = opt$output_directory
)
