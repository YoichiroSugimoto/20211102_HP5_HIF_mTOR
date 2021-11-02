suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
    make_option(
        c("-b", "--bam_file"),
        action="store",
        type = "character",        
        help = "bam file path"),
    make_option(
        c("-p", "--promoter_file"),
        action="store",
        type = "character",        
        help = "promoter location with index (csv file)"),    
    make_option(
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))

splitBamByTssIdx <- function(bam.file, promoters.file, s5.7.7.1.indiv.bams.dir, split.out.bam = 2){
    suppressPackageStartupMessages(library("GenomicAlignments"))
    suppressPackageStartupMessages(library("rtracklayer"))
    suppressPackageStartupMessages(library("data.table"))
    library("magrittr")
    
    promoters.dt <- fread(promoters.file)
    
    promoters.gr <- makeGRangesFromDataFrame(promoters.dt, keep.extra.columns = TRUE)
    
    param <- ScanBamParam(
        what = c("flag", "isize", "seq", "qual")
    )

    bam <- readGAlignmentPairs(bam.file, use.names = TRUE, param = param)

    gal.first <- GenomicAlignments::first(bam)

    gal.tss <- narrow(
        gal.first,
        start = ifelse(strand(gal.first) == "+", 1,  width(gal.first)),
        end = ifelse(strand(gal.first) == "+", 1, width(gal.first))
    )

    sl.idxs <- unique(promoters.dt[, tss_index]) %>% sort

    ol.dt <- data.table(as.data.frame(findOverlaps(gal.tss, promoters.gr)))

    returnAlignedBam <- function(sl.idx, bam, promoters.dt, ol.dt, s5.7.7.1.indiv.bams.dir, split.out.bam){
        sub.idxs <- promoters.dt[tss_index == sl.idx, idx]
        query.idx <- ol.dt[subjectHits %in% sub.idxs, queryHits]
        sl.bam <- bam[query.idx]

        split.out.bam <- ifelse(sl.idx %in% 1:2, split.out.bam, 1) ## Low index does not require splitting
        
        if(length(sl.bam) > 0){

            split.index <- split(1:length(sl.bam), sort(1:length(sl.bam) %% split.out.bam))

            for(i in 1:split.out.bam){            
                out.bam.file <- file.path(
                    s5.7.7.1.indiv.bams.dir,
                    gsub(
                        ".bam$",
                        paste0("_", i, "__idx_", sl.idx, ".bam"),
                        basename(bam.file)
                    )
                )

                export(
                    sl.bam[split.index[[i]]],
                    con = out.bam.file,
                    format = "bam"
                )
            }
        } else {"Nothing to export"}
        
        return()
    }

    temp <- lapply(
        sl.idxs,
        returnAlignedBam,
        bam = bam,
        promoters.dt = promoters.dt,
        ol.dt = ol.dt,
        s5.7.7.1.indiv.bams.dir = s5.7.7.1.indiv.bams.dir,
        split.out.bam = split.out.bam
    )
    return()
}


temp <- splitBamByTssIdx(
    bam.file = opt$bam_file,
    promoters.file = opt$promoter_file,
    s5.7.7.1.indiv.bams.dir = opt$output_directory
)
