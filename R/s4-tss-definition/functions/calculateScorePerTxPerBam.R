suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
    make_option(
        c("-i", "--input_bamfile"),
        action="store",
        type = "character",        
        help = "input bam file"),
    make_option(
        c("-r", "--rdata"),
        action="store",
        type = "character",
        help = "Rdata file for long.scoring.gr and primary.tx.dt"),
    make_option(
        c("-o", "--output_directory"),
        action="store",
        type = "character",
        help = "directory for output files"),
    make_option(
        c("-p", "--filename_postfix"),
        action="store",
        type = "character",
        help = "postfix for output filename")
)

opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("data.table"))

load(opt$rdata)

calculateScorePerTxPerBam <- function(bam.file, long.scoring.gr, all.primary.tx.dt, s4.2.1.2.indiv.tx.score.dir, filename.postfix = "GENCODE"){

    bam <- readGAlignmentPairs(bam.file, use.names = FALSE)

    countOverlapsPerChrs <- function(query, subject){
        ## This funcion is to reduce the memory usage: may not be necessary to use this function with better configuration
        score.dts <- lapply(
            paste0("chr", c(1:22, "X", "Y", "M")),
            function(x, query, subject){
                dt <- data.table(
                    transcript_id = mcols(query[seqnames(query) == x])[, "transcript_id"],
                    count = countOverlaps(
                        query[seqnames(query) == x],
                        subject[seqnames(subject) == x]
                    ),
                    score = mcols(query[seqnames(query) == x])[, "score"]
                )

                score.dt <- dt[, .(sum_score = sum(score * count)), by = transcript_id]
                ## data.table(as.data.frame(
                ##     findOverlaps(
                ##         query[seqnames(query) == x],
                ##         subject[seqnames(subject) == x]
                ##     )
                ## ))
            },
            query = query,
            subject = subject
        )

        score.dt <- rbindlist(score.dts)

        return(score.dt)
    }

    score.dt <- countOverlapsPerChrs(long.scoring.gr, bam)

    ## Check sanity
    if(nrow(score.dt[duplicated(transcript_id)]) != 0){
        stop("Duplicatd tx id existed")
    } else {"OK"}

    tx.score.dt <- merge(
        all.primary.tx.dt[, .(gene_id, gene_name, transcript_id, biotype)],
        score.dt,
        by = "transcript_id"
    )

    tx.score.file <- file.path(
        s4.2.1.2.indiv.tx.score.dir,
        gsub(".bam$", paste0("_", filename.postfix, ".csv"), basename(bam.file))
    )

    fwrite(tx.score.dt, tx.score.file)

    return()

}

temp <- calculateScorePerTxPerBam(
    bam.file = opt$input_bamfile,
    long.scoring.gr = long.scoring.gr,
    all.primary.tx.dt = all.primary.tx.dt,
    s4.2.1.2.indiv.tx.score.dir = opt$output_directory,
    filename.postfix = opt$filename_postfix
)
