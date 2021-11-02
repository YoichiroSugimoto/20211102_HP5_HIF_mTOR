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
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("stringr"))


generateTssBam <- function(sample.name, star.aligned.bam.dir, s2.2.1.tss.bam.dir){
    
    param <- ScanBamParam(
        flag = scanBamFlag(isSecondaryAlignment = FALSE),
        what = c("flag"),
        tag = c("XG") # this tag will define reads containing untemplated G
    )

    input.bam <- list.files(
        star.aligned.bam.dir,
        pattern = glob2rx(paste0(sample.name, "*.bam")),
        full.names = TRUE
    )

    bam <- readGAlignmentPairs(input.bam, use.names = TRUE, param = param)

    gal.first <- GenomicAlignments::first(bam)

    rm(bam)
    gc()
    
    first.bases <- str_split_fixed(names(gal.first), pattern = "\\(|\\)", n = 3)[, 2]
    cigar.first <- cigar(gal.first)

    ## Define reads containing untermplated G at the begining of the reads
    cap.g.flag <- ifelse(
        strand(gal.first) == "+",
        ## strand is plus
                  ifelse(
                      grepl("^1S", cigar.first) & first.bases == "G",
                      1, 0
                  ),
        ## strand is minus
                  ifelse(
                      grepl("1S$", cigar.first) & first.bases == "G",
                      1, 0
                  )    
    )
    mcols(gal.first)$XG <- as.integer(cap.g.flag)

    gal.tss <- c(
        narrow(gal.first[strand(gal.first) == "+"], start = 1, end = 1),
        narrow(gal.first[strand(gal.first) == "-"], start = width(gal.first[strand(gal.first) == "-"]))
    )

    out.tss.bam.file <- file.path(
        s2.2.1.tss.bam.dir,
        paste0(sample.name, ".tss.bam")
    )

    rtracklayer::export(gal.tss, con = out.tss.bam.file, format = "bam")
    
    return()
}


temp <- generateTssBam(
    sample.name = opt$sample_name,
    star.aligned.bam.dir = opt$star_bam_dir,
    s2.2.1.tss.bam.dir = opt$output_directory
)
