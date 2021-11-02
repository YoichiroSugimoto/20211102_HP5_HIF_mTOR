suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(
        c("-i", "--input_bw_file"),
        action="store",
        type = "character",        
        help = "input bam file path"),
    make_option(
        c("-d", "--rdata_for_analysis_range"),
        action="store",
        type = "character",
        help = "Rdata file to define analysis range"),
    make_option(
        c("-o", "--output_dir"),
        action="store",
        type = "character",
        help = "output directory")    
)

opt <- parse_args(OptionParser(option_list=option_list))

createExonCountTable <- function(input.bw.file, analysis.range.rdata, s9.3.1.exon.count.dir){
    
    suppressPackageStartupMessages(require("rtracklayer"))
    suppressPackageStartupMessages(require("magrittr"))
    suppressPackageStartupMessages(require("data.table")) 
   
    load(analysis.range.rdata)
    
    input.bw <- rtracklayer::import(input.bw.file)
    ol.dt <- findOverlaps(
        analysis.range.gr,
        input.bw,
        ignore.strand = TRUE
    ) %>% as.data.frame %>% data.table
    
    ## Data ovelapping with multiple annotations will be ignored
    ol.dt <- ol.dt[!(subjectHits %in% ol.dt[duplicated(subjectHits), subjectHits])]

    ol.sum.dt <- cbind(
        mcols(analysis.range.gr[ol.dt[, queryHits]]) %>% as.data.frame %>% data.table,
        mcols(input.bw[ol.dt[, subjectHits]]) %>% as.data.frame %>% data.table
    )
    ol.sum.count.dt <- ol.sum.dt[, list(score = sum(score)), by = exon_id]
    
    ol.sum.count.dt <- merge(
        data.table(exon_id = mcols(analysis.range.gr)[["exon_id"]]),
        ol.sum.count.dt,
        by = "exon_id",
        all.x = TRUE
    )
    setnafill(ol.sum.count.dt, fill = 0, cols = "score")
    
    setnames(
        ol.sum.count.dt,
        old = "score",
        new = basename(input.bw.file) %>% {gsub("\\.bw", "", .)}
    )
    
    file.path(
        s9.3.1.exon.count.dir,
        paste0(basename(input.bw.file) %>% {gsub("\\.bw", "", .)}, ".csv")
    ) %>%
        {fwrite(ol.sum.count.dt, file = .)}
    
    return()
}

temp <- createExonCountTable(
    input.bw.file = opt$input_bw_file,
    analysis.range.rdata = opt$rdata_for_analysis_range,
    s9.3.1.exon.count.dir = opt$output_dir
)
