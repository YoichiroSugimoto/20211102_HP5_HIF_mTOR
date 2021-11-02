suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(
        c("-i", "--input_bw_file"),
        action="store",
        type = "character",        
        help = "input bam file path"),
    make_option(
        c("-b", "--bedfile_for_analysis_range"),
        action="store",
        type = "character",
        help = "bed file to define analysis range"),
    make_option(
        c("-o", "--output_dir"),
        action="store",
        type = "character",
        help = "output directory")    
)

opt <- parse_args(OptionParser(option_list=option_list))

createDiffBw <- function(input.bw.file, analysis.range.bed, s9.3.1.diff.bw.dir){
    suppressPackageStartupMessages(require("rtracklayer"))
    input.bw <- rtracklayer::import(input.bw.file)
    sl.non.ol.non.first.exon.tx.gr <- rtracklayer::import(analysis.range.bed)
    
    diff.input.bw <- subsetByOverlaps(
        input.bw, sl.non.ol.non.first.exon.tx.gr#, invert = TRUE
    )
    export.bw.file <- file.path(
        s9.3.1.diff.bw.dir,
        basename(input.bw.file)
    )
    rtracklayer::export(diff.input.bw, export.bw.file)
    return()
}


temp <- createDiffBw(
    input.bw.file = opt$input_bw_file,
    analysis.range.bed = opt$bedfile_for_analysis_range,
    s9.3.1.diff.bw.dir = opt$output_dir
)
