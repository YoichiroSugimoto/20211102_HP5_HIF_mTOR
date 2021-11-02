calcMeanCount <- function(base.input.name, count.dt, ref.column.name = "tss_name"){
    require("data.table")
    mean.dt <- count.dt[
       ,
        .(mean_count = rowMeans(.SD)),
        .SDcols = grep(base.input.name, colnames(count.dt), value = TRUE),
        by = ref.column.name
    ]
    
    setnames(mean.dt, old = "mean_count", new = paste0("meanNormCount_", base.input.name))
    return(mean.dt)
}

analyzeDE <- function(total.count.dt,
                      annot.dt,
                      ref.column.name = "tss_name",
                      input.sample.data.df,
                      comparison.name,
                      exp.design,
                      out.dir
                      ){

    require("dplyr")
    require("data.table")
    require("magrittr")
    require("DESeq2")
    
    ## Parse formula
    exp.design.elements  <-
        Reduce(paste, deparse(exp.design)) %>%
        gsub("~", "", .) %>%
        gsub(" ", "", .) %>%
        str_split(., "\\+") %>%
        .[[1]]

    exp.compare.element <- exp.design.elements[length(exp.design.elements)]

    ## Parse comparisons
    base.sample.name <- str_split_fixed(comparison.name, "__", n = 2)[, 1]
    base.comps <- str_split_fixed(comparison.name, "__", n = 2)[, 2] %>%
        str_split_fixed(., "_vs_", n = 2) %>%
        .[1, ]
    base.input.names <- unlist(lapply(
        base.comps,
        function(x){gsub("xx", x, base.sample.name)}
    ))
    
    coef.name <- paste(
        exp.compare.element,
        str_split_fixed(comparison.name, "__", n = 2)[, 2],
        sep = "_"
    )
    
    comparison.regex <- gsub(
        "xx",
        paste0("(", gsub("_vs_", "|", str_split_fixed(comparison.name, "__", n = 2)[, 2] ), ")"),
        str_split_fixed(comparison.name, "__", n = 2)[, 1]
    )

    selected.samples <- rownames(
        input.sample.data.df[grepl(comparison.regex, rownames(input.sample.data.df)), ] 
    )

    dds <- DESeqDataSetFromMatrix(
        countData = deseq2.in.list$count.df[, selected.samples],
        colData = input.sample.data.df[rownames(input.sample.data.df) %in% selected.samples, ],
        design = exp.design
    ) %>%
        DESeq(quiet = TRUE)

    ## Mild pre-filtering is applied
    ## dds <- dds[rowSums(counts(dds)) > 10, ]

    ## CHeck the level of comparisons
    ## print(dds[[comparison.factor]])
    print(paste("Analyzing", comparison.name))
    Reduce(paste, deparse(exp.design)) %>%
        {paste0("With the formula of ", .)} %>%
        print
    print(colData(dds))

    ## Extract no shrinkage fold changes
    raw.fc.res.dt <- DESeq2::results(
                              dds, name = coef.name
                          ) %>%
        as.data.frame %>%
        data.table(keep.rownames = TRUE)
        
    setnames(
        raw.fc.res.dt,
        old = c("rn", "log2FoldChange"),
        new = c(ref.column.name, "log2fc")
    )

    ## Extract shrinkage fold change
    res <- lfcShrink(
        dds,
        type = "apeglm",
        coef = coef.name,
        parallel = TRUE,
        quiet = TRUE
    )

    ## mcols(res, use.names = TRUE)
    plotMA(
        res,
        ylim = c(-2, 2),
        ylab = "log2 fold change",
        main = paste0(comparison.name, "\n mRNA abundance FC in total fraction")    
    )

    res.dt <- data.table(
        as.data.frame(res),
        keep.rownames = TRUE
    ) %>%
        cbind(data.table(dispersion = mcols(dds)[, "dispersion"]))

    setnames(
        res.dt,
        old = c("rn", "log2FoldChange"),
        new = c(ref.column.name, "shrlog2fc")
    )
    
    ## Extract basemean of each sample set
    norm.count.dt <- data.table(counts(dds, normalized = TRUE), keep.rownames = TRUE)
    setnames(norm.count.dt, old = "rn", new = ref.column.name)

    mean.norm.count.dts <- lapply(
        base.input.names,
        calcMeanCount,
        count.dt = norm.count.dt,
        ref.column.name = ref.column.name
    )

    mean.norm.count.dt <- merge(
        mean.norm.count.dts[[1]], mean.norm.count.dts[[2]], by = ref.column.name
    )

    setnames(
        mean.norm.count.dt,
        old = paste0("meanNormCount_", base.input.names),
        new = paste0("meanNormCount_", c("treated", "base"))
    )

    norm.factor <-
        10^6 /
        mean(colSums(mean.norm.count.dt[, .(meanNormCount_treated, meanNormCount_base)]))

    mean.norm.count.dt[, `:=`(
        meanNormCount_treated = meanNormCount_treated * norm.factor,
        meanNormCount_base = meanNormCount_base * norm.factor,
        treated_basename = base.input.names[1],
        base_basename = base.input.names[2]
    )]

    ## Merge all outputs: raw.fc.res.dt, res.dt, mean.norm.count.dt, and annot.dt

    all.res.dt <- merge(
        raw.fc.res.dt[, c(ref.column.name, "log2fc"), with = FALSE],
        res.dt,
        by = ref.column.name
    ) %>%
        merge(
            y = mean.norm.count.dt,
            by = ref.column.name
        ) %>%
        {.[, c(
             ref.column.name,
             "baseMean", "lfcSE", "dispersion", "pvalue", "padj",
             "log2fc", "shrlog2fc",
             "treated_basename", "base_basename",
             "meanNormCount_treated", "meanNormCount_base"
         ), with = FALSE]}
    
    if(is.data.table(annot.dt)){
        all.res.dt <- merge(
            annot.dt,
            all.res.dt,
            by = ref.column.name
        )
    } else {"annotation not added"}
    
    de.table.file <- file.path(
        out.dir,
        paste0(comparison.name, "_DE.csv")
    )

    fwrite(all.res.dt, file = de.table.file)
        
    ## Report summary
    print("Summary of analysis results")
    sum.tb <- table(
        list(
            stat_significance = factor(
                ifelse(
                    all.res.dt[, padj] < 0.1, "padj < 0.1", "not significant"
                ), levels = c("padj < 0.1", "not significant")
            ),
            reg_direction = factor(
                ifelse(all.res.dt[, shrlog2fc] > 0, "Up", "Down"),
                levels = c("Up", "Down")
            )
        )
    )
    print(addmargins(sum.tb))

    return(all.res.dt)
}

print("The following functions are exported: calcMeanCount(), analyzeDE()")
