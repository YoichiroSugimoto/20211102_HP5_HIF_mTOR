suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(
        c("-s", "--sample_name"),
        action="store",
        type = "character",        
        help = "sample name uniquly included in the file name"),    
    make_option(
        c("-i", "--input_directory"),
        action="store",
        type = "character",        
        help = "directory for input files"),
    make_option(
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ShortRead"))

addFirstBaseToFqEntryName <- function(sample.name, processed.fq.step3.dir, processed.fq.step3.2.dir){

    step3.2.in.files <- file.path(
        processed.fq.step3.dir,
        paste0(sample.name, c("_R1.fastq.gz", "_R2.fastq.gz"))
    )
    
    step3.2.out.files <- file.path(
        processed.fq.step3.2.dir,
        basename(step3.2.in.files)
    )

    ## Process first read
    fq1 <- readFastq(step3.2.in.files[1])

    id.chars <- as.character(id(fq1))
    splitted.id.chars <- str_split_fixed(id.chars, pattern = "_", n = 2)

    first.bases <- as.character(narrow(sread(fq1), start = 1, end = 1))

    new.id.chars <- paste0(
        splitted.id.chars[, 1],
        ":(",
        first.bases,
        ")_",
        splitted.id.chars[, 2]
    )

    fq1 <- ShortReadQ(sread = sread(fq1), id = BStringSet(new.id.chars), quality = quality(fq1))

    if(file.exists(step3.2.out.files[1])){
        file.remove(step3.2.out.files[1])
    } else {}

    writeFastq(fq1, file = step3.2.out.files[1])

    ## Process second read
    fq2 <- readFastq(step3.2.in.files[2])
    fq2 <- ShortReadQ(sread = sread(fq2), id = BStringSet(new.id.chars), quality = quality(fq2))
    
    if(file.exists(step3.2.out.files[2])){
        file.remove(step3.2.out.files[2])
    } else {}

    writeFastq(fq2, file = step3.2.out.files[2])
    
    return()

}

addFirstBaseToFqEntryName(
    sample.name = opt$sample_name,
    processed.fq.step3.dir = opt$input_directory,
    processed.fq.step3.2.dir = opt$output_directory
)
