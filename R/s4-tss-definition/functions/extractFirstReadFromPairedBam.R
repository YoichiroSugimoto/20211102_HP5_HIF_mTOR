## This function trim soft clipped bases because of the reverse transcription of mRNA cap from bam files

suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
    make_option(
        c("-i", "--input_bamfile"),
        action="store",
        type = "character",        
        help = "input bam file"),
    make_option(
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("rtracklayer"))

extractFirstReadFromPairedBam <- function(bam.file, s5.1.4.b.merged.first.strand.bam.dir){

    what <- c("flag", "seq", "qual") # Optional column not necessary?
    param <- ScanBamParam(what = what, tag = c("NH", "HI", "AS", "nM"))
    bam <- readGAlignmentPairs(bam.file, use.names = TRUE, param = param)

    gal.first <- GenomicAlignments::first(bam)
    rm(bam)
    gc()

    names(gal.first) <- str_split_fixed(names(gal.first), "_", n = 2)[, 1]

    out.bam.file <- file.path(
        s5.1.4.b.merged.first.strand.bam.dir,
        gsub(".bam$", "-first-s.bam", basename(bam.file))
    )

    export(gal.first, con = out.bam.file, format = "bam")

    return()
}

temp <- extractFirstReadFromPairedBam(opt$input_bamfile, opt$output_directory)

