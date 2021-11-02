## ## The following data are for development
## library("stringr")
## library("dplyr")
## library("data.table")
## results.dir <- file.path("../../../results")
## ## results.dir <- file.path("../../../../20200517_five_prime_seq_for_VHL_loss/results")
## s8.dir <- file.path(results.dir, "s8-analysis-of-translation")
## s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs")
## s8.1.rdata <- file.path(s8.1.dir, "s8.1.test.analyzeDtg.Rdata")
## load(s8.1.rdata)
## s8.1.dir <- file.path(s8.dir, "s8-1-differentially-translated-mRNAs", "gene-level-dte")
## count.df = deseq2.poly.in.list$count.df
## annot.df = deseq2.poly.in.list$annot.df
## ref.column.name = "gene_id"
## sizefactor.col = "size_factor"
## comparison.name = translation.comparison.dt[3, comparison]
## deseq2.formula = translation.comparison.dt[, exp_formula][[3]]
## input.sample.data.df = subsetColdata(comparison.name, poly.coldata.df)
## s3.alignment.stats.dir <- file.path(results.dir, "s3-alignment-statistics")
## s3.4.poly.size.factor.dir <- file.path(s3.alignment.stats.dir, "polysome_size_factor")
## poly.sizefactor.dt <- fread(
##     file = file.path(
##         s3.4.poly.size.factor.dir,
##         "library_size_factor_by_ERCC.csv"
##     )
## )

subsetColdata <- function(comparison.name, poly.coldata.df){

    regex.comp <- str_split_fixed(comparison.name, "__", n = 2)[, 2] %>%
        gsub("_vs_", "$|^", .) %>%
        paste0("(^", ., "$)")
    
    regexElement <- function(element, regex.comp){
        element <- case_when(
            element == "xx" ~ regex.comp,
            element == "yy" ~ "*",
            TRUE ~ paste0("^", element, "$")
        )
        return(element)
    }

    cell <- str_split_fixed(comparison.name, pattern = "_", n = 6)[, 1] %>%
        regexElement(., regex.comp)
    VHL <- str_split_fixed(comparison.name, pattern = "_", n = 6)[, 2] %>%
        regexElement(., regex.comp)
    EIF4E2 <- str_split_fixed(comparison.name, pattern = "_", n = 6)[, 3] %>%
        regexElement(., regex.comp)
    clone <- str_split_fixed(comparison.name, pattern = "_", n = 6)[, 4] %>%
        regexElement(., regex.comp)
    treatment <- str_split_fixed(comparison.name, pattern = "_", n = 6)[, 5] %>%
        regexElement(., regex.comp)
    
    rn.dt <- str_split_fixed(rownames(poly.coldata.df), "_", n = 8) %>%
        data.table %>%
        cbind(
            data.table(rn = rownames(poly.coldata.df)),
            .
        )

    sl.rownames <- rn.dt[
        grepl(cell, V1) &
        grepl(VHL, V2) &
        grepl(EIF4E2, V3) &
        grepl(clone, V5) &
        grepl(treatment, V6),
        rn
    ]

    return(poly.coldata.df[rownames(poly.coldata.df) %in% sl.rownames, ])
}

## Define functions for MRL calculation
normalizeHP5readsByStandard <- function(polysome.count.df, polysome.sizefactor.dt, sizefactor.col, ref.column.name){

    if(all(colnames(polysome.count.df) != polysome.sizefactor.dt[, sample_name])){
        stop("Samples in the count table and size factor table do not match")
    } else {"OK"}

    int.normalized.polysome.count.df <- sweep(
        x = polysome.count.df,
        MARGIN = 2,
        STATS = polysome.sizefactor.dt[, get(sizefactor.col)],
        FUN = "/"
    )

    ## Confirm the normalization
    if(!isTRUE(all.equal(
            target = colSums(polysome.count.df) / colSums(int.normalized.polysome.count.df),
            current = polysome.sizefactor.dt[, get(sizefactor.col)],
            tolerance = 1.0e-5,
            check.attributes = FALSE
        ))){
        stop("Check normalization")
    } else {"OK"}

    ## Normalize so that the average count per fraction becomes million
    normalized.polysome.count.df <- sweep(
        x = int.normalized.polysome.count.df,
        MARGIN = 2,
        STATS = colSums(int.normalized.polysome.count.df) %>% mean %>% {10^6/.},
        FUN = "*"
    )

    print("QC stats for library size normalization")
    data.table(
        polysome.sizefactor.dt[, c("sample_name" , sizefactor.col), with = FALSE],
        original_libsize = colSums(polysome.count.df),
        int_libsize = colSums(polysome.count.df) / polysome.sizefactor.dt[, get(sizefactor.col)],
        norm_libsize = colSums(normalized.polysome.count.df)
    ) %>%
        print()

    normalized.polysome.count.dt <- data.table(
        normalized.polysome.count.df,
        keep.rownames = TRUE
    )
    setnames(normalized.polysome.count.dt, old = "rn", new = ref.column.name)

    return(normalized.polysome.count.dt)
}


calculateWeightedMRL <- function(normalized.polysome.count.dt, sl.ribo.sample.groups, ref.column.name){
    
    calculateWeightedMRLperGroup <- function(ribo.sample.group, count.dt){
        sl.count.dt <- count.dt[
          , c(
                ref.column.name,
                grep(ribo.sample.group, colnames(count.dt), value = TRUE)),
            with = FALSE]

        sl.count.dt[, `:=`(
            weightedRLR = rowSums(sweep(
                .SD,
                2,
                str_extract(
                    colnames(sl.count.dt)[grepl(ribo.sample.group, colnames(sl.count.dt))],
                    "ribo[[:digit:]]"
                ) %>%
                {gsub("ribo", "", .)} %>%
                as.integer
               ,
                FUN = "*"
            )),
            totalRLR = rowSums(.SD)
        ),
        .SDcols = grep(ribo.sample.group, colnames(sl.count.dt), value = TRUE)
        ]

        data.cols <- c("weightedRLR", "totalRLR")
        
        sl.count.dt <- sl.count.dt[, c(ref.column.name, data.cols), with = FALSE]
        setnames(
            sl.count.dt,
            old = data.cols,
            new = paste0(gsub("_ribo\\[\\[:digit:\\]\\]\\(\\|A\\|B\\)", "", ribo.sample.group), "_", data.cols)
        )
        setkeyv(sl.count.dt, ref.column.name)
        return(sl.count.dt)
    }
    
    mrl.dts <- lapply(
        sl.ribo.sample.groups,
        calculateWeightedMRLperGroup,
        count.dt = normalized.polysome.count.dt
    )

    mrl.count.dt <- Reduce(
        function(...) merge(..., all = TRUE, by = ref.column.name),
        mrl.dts
    )

    return(mrl.count.dt)
}

calcMrlStat <- function(mrl.count.dt, annot.df, ref.column.name, base.input.names){

    mrl.summary.dt <- copy(mrl.count.dt)

    m.mrl.summary.dt <- melt(
        mrl.summary.dt,
        id.vars = ref.column.name,
        value.name = "TorW_RLR",
        variable.name = "full_sample_name"
    ) %>%
        {.[, `:=`(
             sample_name = gsub("(_weightedRLR$|_totalRLR$)", "", full_sample_name),
             RLR_type = str_extract(full_sample_name, "(weightedRLR$|totalRLR$)")
         )]}

    m.mrl.summary.dt[, `:=`(
        input_basename = gsub("EIF4E2_[[:alnum:]]*_", "EIF4E2_", sample_name) %>%
            {case_when(
                 grepl(
                     gsub("_yy", "_*[[:alnum:]]", base.input.names[1]),
                     .
                 ) ~ base.input.names[1],
                 grepl(
                     gsub("_yy", "_*[[:alnum:]]", base.input.names[2]),
                     .
                 ) ~ base.input.names[2]
             )} %>%
            factor(levels = base.input.names)
    )]

    ## Sanity check
    if(nrow(m.mrl.summary.dt[is.na(input_basename)]) != 0){
        stop("base_inputname does not match with some of sample names")
    } else {"OK"}
    
    d.mrl.summary.dt <- dcast(
        m.mrl.summary.dt,
        formula = as.formula(
            paste(ref.column.name, "+ sample_name + input_basename ~ RLR_type")
        ),
        value.var = "TorW_RLR"
    )

    d.mrl.summary.dt[, raw_MRL := weightedRLR / totalRLR]

    mrl.stat.dt <- d.mrl.summary.dt[, list(
        MRL = mean(raw_MRL),
        MRLsd = sd(raw_MRL),
        totalRLR = mean(totalRLR),
        totalRLRsd = sd(totalRLR),
        weightedRLR = mean(weightedRLR),
        weightedRLRsd = sd(weightedRLR)
    ), by = list(input_basename, get(ref.column.name))]
    
    setnames(mrl.stat.dt, old = "get", new = ref.column.name)
    
    mrl.stat.dt <- mrl.stat.dt %>%
        {.[!duplicated(paste(get(ref.column.name), input_basename))]} %>%
        dcast(
            formula = as.formula(paste(ref.column.name,  "~ input_basename")),
            value.var = c(
                "MRL", "MRLsd",
                "totalRLR", "totalRLRsd",
                "weightedRLR", "weightedRLRsd"
            )
        ) %>%
        {.[, `:=`(
             treated_basename = base.input.names[1],
             base_basename = base.input.names[2]
         )]}

    setnames(
        mrl.stat.dt,
        old = colnames(mrl.stat.dt),
        new = colnames(mrl.stat.dt) %>%
            {gsub(base.input.names[1], "treated", .)} %>%
            {gsub(base.input.names[2], "base", .)}
    )

    mrl.stat.dt[, `:=`(
        MRL_log2fc = log2(MRL_treated / MRL_base),
        totalRLR_log2fc = log2(totalRLR_treated / totalRLR_base),
        weightedRLR_log2fc = log2(weightedRLR_treated / weightedRLR_base)
    )]
    
    mrl.count4export.dt <- merge(
        data.table(annot.df),
        mrl.stat.dt,
        by = ref.column.name,
        all = TRUE
    ) ## %>%
    ##   merge(
    ##       y = mrl.count.dt,
    ##       by = ref.column.name,
    ##       all = TRUE
    ##   )
    return(mrl.count4export.dt)
    
}

analyzeDtg <- function(
                       count.df,
                       annot.df,
                       ref.column.name = "tss_name",
                       comparison.name,
                       deseq2.formula,
                       input.sample.data.df,
                       poly.sizefactor.dt,
                       sizefactor.col = "size_factor",
                       s8.1.dir,
                       processors
                       ){

    library("tidymodels")
    library("DESeq2")
    library("dplyr")
    library("stringr")
    library("data.table")
    library("magrittr")
    library("matrixStats")
    library("ggplot2")
    library("pheatmap")
    
    test.fractions <- paste0("ribo", 1:8)

    ## Prepare reduced formula
    deseq2.reduced.formula.elements <- unlist(str_split(as.character(deseq2.formula)[2], " \\+ "))

    deseq2.reduced.formula <- as.formula(
        paste("~", paste(
                       deseq2.reduced.formula.elements[
                           1:(length(deseq2.reduced.formula.elements) - 1)
                       ],
                       collapse = " + ")
              )
    )

    ## Extract compared factors
    base.sample.name <- str_split_fixed(comparison.name, "__", n = 2)[, 1]
    base.comps <- str_split_fixed(comparison.name, "__", n = 2)[, 2] %>%
        str_split_fixed(., "_vs_", n = 2) %>%
        .[1, ]
    base.input.names <- unlist(lapply(
        base.comps,
        function(x){gsub("xx", x, base.sample.name)}
    ))
    
    ## Prepare factor to compare
    compared.factor <- deseq2.reduced.formula.elements[
        length(deseq2.reduced.formula.elements)
    ] %>%
        gsub("fraction:", "", .)

    ## Prepare coefccient for DESeq2
    comparison.elements <- str_split_fixed(comparison.name, "__", n = 2)[, 2] %>%
        str_split_fixed(., "_vs_", n = 2)
    deseq2.coef <- paste0(".", compared.factor, comparison.elements[, 1])
    
    print("The data that will be used for this analysis")
    print(input.sample.data.df)

    ## Comparison name for export
    comparison.name4export <- comparison.name %>%
        gsub("\\(", "", .) %>%
        gsub("\\)", "", .) %>%
        gsub("\\|", "-", .)

    ## subset count.df
    test.count.df <- count.df[, rownames(input.sample.data.df)]

    ## subset size.factor.dt
    poly.sizefactor.dt <- poly.sizefactor.dt[
       , sample_name := gsub("polysome_", "", sample_name)
        ] %>%
        {.[
            sample_name %in% rownames(input.sample.data.df)
        ]}
        
    ## #######################################################################
    ## Mean ribosome load calcuation
    print("MRL calculation")

    ## Here I only analyze polysome fractionated samples
    polysome.fractionated.samples <- input.sample.data.df[
        input.sample.data.df$fraction %in% test.fractions,
        ] %>% rownames

    polysome.count.df <- test.count.df[
       , colnames(test.count.df) %in% polysome.fractionated.samples
    ]

    polysome.sizefactor.dt <- poly.sizefactor.dt[
        sample_name %in% polysome.fractionated.samples
    ][
        order(match(sample_name, polysome.fractionated.samples))
    ]

    normalized.polysome.count.dt <- normalizeHP5readsByStandard(
        polysome.count.df = polysome.count.df,
        polysome.sizefactor.dt = polysome.sizefactor.dt,
        sizefactor.col = sizefactor.col,
        ref.column.name = ref.column.name
    )
        
    sl.ribo.sample.groups <- grep(
        "ribo1", rownames(input.sample.data.df), value = TRUE) %>%
        gsub("ribo1", "ribo[[:digit:]](|A|B)", .)


    mrl.count.dt <- calculateWeightedMRL(
        normalized.polysome.count.dt = normalized.polysome.count.dt,
        sl.ribo.sample.groups = sl.ribo.sample.groups,
        ref.column.name = ref.column.name
    )
        
    mrl.count4export.dt <- calcMrlStat(
        mrl.count.dt = mrl.count.dt,
        annot.df = annot.df,
        ref.column.name = ref.column.name,
        base.input.names = base.input.names
    )
    
    mrl.for.export.file <- file.path(
        s8.1.dir,
        paste0(comparison.name4export, "-mean_ribosome_loading.csv")
    )
    
    fwrite(mrl.count4export.dt, file = mrl.for.export.file)

    ############################################################
    ## DESeq2 interaction model: translation regulation classification
    print("LRT comparison")

    Reduce(paste, deparse(deseq2.formula)) %>%
        paste0("Formula for LRT comparison: ", .) %>%
        print
    Reduce(paste, deparse(deseq2.reduced.formula)) %>%
        paste0("Reduced formula for LRT comparison: ", .) %>%
        print

    lrt.test.count.df <- test.count.df[
      , str_split_fixed(colnames(test.count.df), "_", n = 8)[, 7] %in% c(test.fractions)
    ]
        
    lrt.input.sample.data.df <- input.sample.data.df[colnames(lrt.test.count.df), ]
    
    lrt.dds <- DESeqDataSetFromMatrix(
        countData = lrt.test.count.df,
        colData = lrt.input.sample.data.df,
        design = deseq2.formula
    ) %>%
        estimateSizeFactors
    
    lrt.dds <- DESeq(
        lrt.dds,
        test = "LRT",
        reduced = deseq2.reduced.formula,
        quiet = TRUE
    )
    
    lrt.res.dt <- results(lrt.dds) %>% as.data.frame %>%
        data.table(keep.rownames = ref.column.name)

    setnames(
        lrt.res.dt,
        old = "padj",
        new = "padj_translation"
    )

    dtg.dt <-  merge(
        lrt.res.dt[, c(ref.column.name, "padj_translation"), with = FALSE],
        mrl.count4export.dt[
          , c(ref.column.name, "MRL_log2fc", "MRL_treated", "MRL_base"), with = FALSE
        ],
        by = ref.column.name
    )

    sl.dtg.dt <- dtg.dt[padj_translation < 0.1 & !is.na(MRL_log2fc)]

    if(nrow(sl.dtg.dt) > 0){

        ## kmean
        ## sl.dtg.dt[, kmeans_cluster := kmeans(MRL_log2fc,centers = 2)$cluster - 1]
        
        
        ## is.zero.higher.mrl <-
        ##     mean(sl.dtg.dt[kmeans_cluster == 0, MRL_log2fc]) >
        ##     mean(sl.dtg.dt[kmeans_cluster == 1, MRL_log2fc])

        ## sl.dtg.dt[
        ##   , translational_regulation := case_when(
        ##         is.zero.higher.mrl & kmeans_cluster == 0 ~ "Up",
        ##         is.zero.higher.mrl & kmeans_cluster == 1 ~ "Down",
        ##         !is.zero.higher.mrl & kmeans_cluster == 0 ~ "Down",
        ##         !is.zero.higher.mrl & kmeans_cluster == 1 ~ "Up"
        ##     )
        ## ]

        ## Gausian mixture models
        ## library("mixtools")
        ## mix.model <- normalmixEM(sl.dtg.dt[, MRL_log2fc], k = 2)
        
        ## sl.dtg.dt <- cbind(
        ##     sl.dtg.dt, data.table(mix.model$posterior)
        ## )

        ## setnames(
        ##     sl.dtg.dt,
        ##     old = c("comp.1", "comp.2"),
        ##     new = sort(c("p_up", "p_down"), decreasing = mix.model$mu %>% {.[1] < .[2]})
        ## )
        
        ## sl.dtg.dt[
        ##    ,
        ##     translational_regulation := case_when(
        ##         p_down > 0.9 ~ "Down",
        ##         p_up > 0.9 ~ "Up",
        ##         TRUE ~ "Unclassified"
        ##     )
        ## ]

        ## Very simple
        sl.dtg.dt[
          , translational_regulation := case_when(
                MRL_log2fc < dtg.dt[, median(MRL_log2fc, na.rm = TRUE)] ~ "Down",
                MRL_log2fc > dtg.dt[, median(MRL_log2fc, na.rm = TRUE)] ~ "Up",
                TRUE ~ "Unclassified"
            )
        ]

        dtg.dt <- merge(
            dtg.dt,
            sl.dtg.dt[, c(ref.column.name, "translational_regulation"), with = FALSE],
            by = ref.column.name,
            all = TRUE
        )

        dtg.dt[, translational_regulation := case_when(
                     is.na(translational_regulation) ~ "Unclassified",
                     TRUE ~ translational_regulation
                 )]

        g1 <- ggplot(
            data = dtg.dt[!is.na(padj_translation)],
            aes(
                x = MRL_log2fc,
                fill = translational_regulation
            )
        ) +
            geom_histogram(binwidth = 0.025, alpha = 0.5, position = "identity") +
            geom_vline(
                xintercept = dtg.dt[
                    !is.na(padj_translation), median(MRL_log2fc, na.rm = TRUE)
                ]
            ) +
            facet_grid(
                ifelse(padj_translation < 0.1, "Sig", "Non sig") ~ .,
                scales = "free"
            ) +
            ggsci::scale_fill_aaas()

        print(g1)
        
    } else {
        dtg.dt[, translational_regulation := "Unclassified"]
    }
    
    ## Export normalized count
    tpm.poly.count4export.dt <- merge(
        data.table(annot.df),
        normalized.polysome.count.dt,
        by = ref.column.name,
        all = TRUE
    )
    
    count.for.export.file <- file.path(
        s8.1.dir,
        paste0(comparison.name4export, "-normalized_count.csv")
    )

    fwrite(tpm.poly.count4export.dt, file = count.for.export.file)

    ## Export mean ribosome loading
    ## 1st, calculate MRL and the sd
    
    ## Export all results
    dtg4export.dt <- merge(
        data.table(annot.df),
        dtg.dt,
        by = ref.column.name,
        all = TRUE
    )
    
    fwrite(dtg4export.dt, file = file.path(s8.1.dir, paste0(comparison.name4export, ".csv")))

    return()
}

print("The following functions were exported: analyzeDtg(), subsetColdata()")
