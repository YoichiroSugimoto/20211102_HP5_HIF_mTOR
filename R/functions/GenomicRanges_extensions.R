## The functions here are to extend GenomicRanges

countPerRangeFromBams <- function(input.bam.files, gr, out.file.name, processors, range.name = "tss_name", bam.file.prefix = ".tss.bam$", sample.names = NULL){
    ## This function count the number of TSS bam mapped within ranges specified by genomic ranges
    ## input.bam.files: vector of input bam files
    ## ge: GRanges of range
    ## out.file.name: output count table filename
    ## processors: the number of processors to be used
    ## range.name: the names of range
    ## bam.file.prefix: input bamfile is expected to have the name of sample.name + bam.file.prefix
    ## sample.names: if sample.names are specified here, the order of column will be sorted in this order
    require("GenomicAlignments")
    
    createCountTable <- function(bam.file, bam.file.prefix){
        bam <- readGAlignments(bam.file)
        ol.count <- countOverlaps(gr, bam)
        ol.count.dt <- data.table(V1 = mcols(gr)[, range.name], count = ol.count)
        ol.count.dt <- ol.count.dt[, .(count = sum(count)), by = V1]
        setnames(
            ol.count.dt,
            old = c("V1", "count"),
            new = c(range.name, gsub(bam.file.prefix, "", basename(bam.file)))
        )
        setkeyv(ol.count.dt, range.name)
        return(ol.count.dt)
    }
    
    count.dts <- mclapply(
        input.bam.files,
        createCountTable,
        bam.file.prefix = bam.file.prefix,
        mc.cores = processors
    )
    
    count.dt <- Reduce(function(...) merge(..., by = range.name), count.dts)
    setkeyv(count.dt, range.name)

    if(!is.null(sample.names)){ 
        count.dt <- count.dt[, c(range.name, sample.names), with = FALSE]
    } else {"No columnname sorting"}
    
    fwrite(count.dt, file = out.file.name)
    
    return(count.dt)
}
