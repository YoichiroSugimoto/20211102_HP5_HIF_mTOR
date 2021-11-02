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

trimSoftClippedGFromBam <- function(bam.file, leaf.no.softclip.bam.dir){

    what <- c("flag", "seq") # Optional column not necessary?
    param <- ScanBamParam(what = what, tag = c("NH", "HI", "AS", "nM"))
    bam <- readGAlignmentPairs(bam.file, use.names = TRUE, param=param)

    gal.first <- GenomicAlignments::first(bam)
    gal.second <- GenomicAlignments::last(bam, real.strand = FALSE)

    rm(bam)

    cigar.first <- cigar(gal.first)

    cap.g.flag <- ifelse(
        strand(gal.first) == "+",
                  ifelse(
                      grepl("^1S", cigar.first) &
                      substr(elementMetadata(gal.first)$seq, start = 1, stop = 1) == "G",
                      TRUE, FALSE),
                  ifelse(
                      grepl("1S$", cigar.first) &
                      substr(elementMetadata(gal.first)$seq,
                             start = nchar(elementMetadata(gal.first)$seq),
                             stop = nchar(elementMetadata(gal.first)$seq)) == "C",
                      TRUE, FALSE)    
    )

    trimmed.gal.first <- qnarrow(
        gal.first,
        start = ifelse(cap.g.flag, 2, 1)
    )

    mcols(trimmed.gal.first)$seq <- NULL
    mcols(gal.second)$seq <- NULL
    
    trimmed.galp <- GAlignmentPairs(
        trimmed.gal.first,
        gal.second,
        names = names(trimmed.gal.first)
    )

    output.bam <- file.path(
        leaf.no.softclip.bam.dir,
        basename(bam.file)
    )

    export(trimmed.galp, con = output.bam, format = "bam")

    rm(gal.first, gal.second, trimmed.galp)
    gc()
    
    return()
}

temp <- trimSoftClippedGFromBam(opt$input_bamfile, opt$output_directory)
