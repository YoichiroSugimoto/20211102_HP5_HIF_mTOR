
## This function filter TSS with expressing less than `threshold.ratio` to the most highly expressed TSS in a gene
thresholdTssByRatio <- function(tss.dt, threshold.ratio = 0.1, tss.name = "tss_name", ref.expression.col = "baseMean_cyto"){
    ## tss.dt <- data.table(gene_id = c("a", "a", "a", "a"), baseMean_cyto = c(10, 5, 0.1, 0.1), baseMean_nuc = c(10, 0.1, 5, 0.1), tss_name = paste0("TSS_", 1:4))
    require("data.table")
    temp.tss.dt <- copy(tss.dt)
    temp.tss.dt[,  gene_max := max(get(ref.expression.col)), by = gene_id]
    temp.tss.dt <- temp.tss.dt[
        get(ref.expression.col) > threshold.ratio * gene_max
    ]
    tss.dt <- tss.dt[get(tss.name) %in% temp.tss.dt[, get(tss.name)]]
    return(tss.dt)
}

## Count bt tss_name to count by gene_id
countByGeneFromTss <- function(count.dt){
    count.dt[, gene_id := str_split_fixed(tss_name, "_", n = 2)[, 1]]
    count.per.gene.dt <- count.dt[
        grepl("^ENSG", gene_id), lapply(.SD, sum), by = gene_id,
        .SDcols = colnames(count.dt)[!(colnames(count.dt) %in% c("tss_name", "gene_id"))]
    ]

    return(count.per.gene.dt)
}

## This function filters polysome data
filterPolysomeByCount <- function(match.pattern, count.dt, ref.colname = "tss_name", min.count.th = 10, min.datapoint.th = 6){

    count.dt <- copy(count.dt[, c(
        ref.colname,
        grep(match.pattern, colnames(count.dt), value = TRUE)
    ), with = FALSE])

    nonzero.col.num <- count.dt[
      , grep(match.pattern, colnames(count.dt), value = TRUE), with = FALSE
    ] %>%
        {rowSums(. > min.count.th)}
    ## Conditions for VHL data will be used only when necessary
    filtered.gene.id <- count.dt[nonzero.col.num > min.datapoint.th, get(ref.colname)]
    print(
        paste0(
            "The number of genes/transcripts after minimum count across fraction filtering: ",
            length(filtered.gene.id))
    )
    return(filtered.gene.id)
}
