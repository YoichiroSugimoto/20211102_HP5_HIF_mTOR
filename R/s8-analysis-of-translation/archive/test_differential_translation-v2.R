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
## comparison.name = translation.comparison.dt[1, comparison]
## deseq2.formula = translation.comparison.dt[, exp_formula][[1]]
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

    ## Calculate log2 fold change of each fraction
    lrt.vsd <- vst(lrt.dds, blind = TRUE)
    lrt.vsd.dt <- assay(lrt.vsd) %>% data.table(keep.rownames = ref.column.name)

    m.lrt.vsd.dt <- melt(
        lrt.vsd.dt,
        id.vars = ref.column.name,
        value.name = "VST",
        variable.name = "sample_name"
    ) 

    m.lrt.vsd.dt[, `:=`(
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
            factor(levels = base.input.names),
        fraction = paste0(
            "log2fc_",
            str_split_fixed(sample_name, "_", n = 8)[, 7]
            )
    )]

    ## Sanity check
    if(nrow(m.lrt.vsd.dt[is.na(input_basename)]) != 0){
        stop("base_inputname does not match with some of sample names")
    } else {"OK"}
    
    m.lrt.vsd.stat.dt <- m.lrt.vsd.dt[, list(
        VST_mean = mean(VST),
        VST_sd = sd(VST)
    ), by = list(input_basename, get(ref.column.name), fraction)]

    setnames(m.lrt.vsd.stat.dt, old = "get", new = ref.column.name)

    d.lrt.vsd.stat.dt <- dcast(
        m.lrt.vsd.stat.dt,
        paste0(ref.column.name, " + fraction ~ input_basename"),
        value.var = "VST_mean"
    )

    d.lrt.vsd.stat.dt[, fraction_log2fc := 
                            get(base.input.names[1]) - get(base.input.names[2])
                      ]

    dd.lrt.vsd.stat.dt <- dcast(
        d.lrt.vsd.stat.dt,
        paste0(ref.column.name, " ~ fraction"),
        value.var = "fraction_log2fc"
    )
    
    dtg.dt <- Reduce(
        function(...) merge(..., all = TRUE, by = ref.column.name),
        c(
            list(
                lrt.res.dt[, c(ref.column.name, "padj_translation"), with = FALSE],
                dd.lrt.vsd.stat.dt
            ),
            list(mrl.count4export.dt[, c(ref.column.name, "MRL_log2fc", "MRL_treated", "MRL_base"), with = FALSE])
        )
    )

    dtg.dt <- dtg.dt[, c(
        ref.column.name,
        "padj_translation", paste0("log2fc_", test.fractions),
        c("MRL_log2fc", "MRL_treated", "MRL_base")
    ), with = FALSE]
    
    if(nrow(dtg.dt[padj_translation < 0.1]) > 10){
        ## First classify translational regulation direction using kmeans and logistic regression
        dte.classification.dt <- copy(dtg.dt[!is.na(MRL_log2fc)])[
            padj_translation < 0.1,
            c(ref.column.name, paste0("log2fc_ribo", 1:8), "MRL_log2fc"),
            with = FALSE
        ] %>%
            {.[, kmeans_cluster := kmeans(
                     .[, c("MRL_log2fc"), with = FALSE],
                     centers = 2
                 )$cluster - 1]}

        is.zero.higher.mrl <-
            mean(dte.classification.dt[kmeans_cluster == 0, MRL_log2fc]) >
            mean(dte.classification.dt[kmeans_cluster == 1, MRL_log2fc])
        
        ## Rename kmeans cluster so that translation up-regulated group is 1
        dte.classification.dt[
          , kmeans_cluster := case_when(
                is.zero.higher.mrl & kmeans_cluster == 0 ~ "Up",
                is.zero.higher.mrl & kmeans_cluster == 1 ~ "Down",
                !is.zero.higher.mrl & kmeans_cluster == 0 ~ "Down",
                !is.zero.higher.mrl & kmeans_cluster == 1 ~ "Up"
            ) %>% factor(levels = c("Down", "Up"), ordered = TRUE)]

        g.lr.1 <- ggplot(
            data = dte.classification.dt,
            aes(
                x = MRL_log2fc,
                fill = kmeans_cluster
            )
        ) +
            geom_histogram(
                alpha = 0.6, binwidth = 0.05#,
                ## position = "identity"
            ) +
            khroma::scale_fill_bright() +
            ggtitle("Translational regulation grouping by kmeans") +
            theme(
                legend.position = "bottom"
            )
        print(g.lr.1)
        
        logfit.translation <- logistic_reg() %>% 
            set_engine("glm") %>% 
            set_mode("classification") %>% 
            parsnip::translate() %>%
            fit(
                kmeans_cluster ~
                    log2fc_ribo1 + log2fc_ribo2 + log2fc_ribo3 + log2fc_ribo4 +
                    log2fc_ribo5 + log2fc_ribo6 + log2fc_ribo7 + log2fc_ribo8,
                data = as.data.frame(dte.classification.dt)
            )
        print(tidy(logfit.translation))

        dte.classification.dt <- dte.classification.dt %>%
            cbind(
                predict(
                    logfit.translation, new_data = dte.classification.dt,
                    type = "prob"
                ),
                predict(
                    logfit.translation, new_data = dte.classification.dt
                )
            ) %>%
            {.[, `:=`(
                 prediction_prob = pmax(get(".pred_Down"), get(".pred_Up")),
                 prediction_p = 1 - pmax(get(".pred_Down"), get(".pred_Up"))
             )]} %>%
            {.[, prediction_pseudoFDR := p.adjust(prediction_p, method = "fdr")]}
        
        print("metics for logistic regression")
        dte.classification.dt %>%
            conf_mat(truth = kmeans_cluster, estimate = .pred_class) %>% print
        dte.classification.dt %>%
            metrics(truth = kmeans_cluster, estimate = .pred_class) %>% print

        dte.classification.dt[
          , translational_regulation := case_when(
                get(".pred_class") == "Up" & prediction_prob > 0.9 ~ "Up",
                get(".pred_class") == "Down" & prediction_prob > 0.9 ~ "Down",
                TRUE ~ "Uncertain"
            )]
        
        print("After requiring 90% probability by logistic regression")
        table(dte.classification.dt[, translational_regulation]) %>% print
        
        g.lr.2 <- ggplot(
            data = dte.classification.dt,
            aes(
                x = MRL_log2fc,
                fill = translational_regulation
            )
        ) +
            geom_histogram(alpha = 0.6, position = "identity", binwidth = 0.05) +
            khroma::scale_fill_bright() +
            theme(
                legend.position = "bottom"
            ) +
            ggtitle("After logistic regression classification")
        print(g.lr.2)
        
        fc.mat <- as.matrix(
            dtg.dt[padj_translation < 0.1, paste0("log2fc_", test.fractions), with = FALSE],
            rownames = unlist(dtg.dt[padj_translation < 0.1, ref.column.name, with = FALSE])
        )

        cal_z_score <- function(x){
            (x - mean(x)) / sd(x)
        }

        fc.mat <- fc.mat[complete.cases(fc.mat), ]
        scaled.fc.mat <- apply(fc.mat, 1, cal_z_score)

        c <- cor(scaled.fc.mat, method = "spearman") 
        d <- as.dist(1-c)
        
        pheatmap.annot.df <- dte.classification.dt[
            order(match(get(ref.column.name), colnames(scaled.fc.mat))),
            translational_regulation
        ] %>% as.data.frame
        rownames(pheatmap.annot.df) <- dte.classification.dt[
            order(match(get(ref.column.name), colnames(scaled.fc.mat))),
            get(ref.column.name)
        ]
        
        pheatmap(
            t(scaled.fc.mat),
            cluster_cols = FALSE,
            clustering_distance_rows = d,
            clustering_method = "ward.D2",
            annotation_row = pheatmap.annot.df,
            border_color = FALSE,
            show_rownames = FALSE
        )    

        dtg.dt <- merge(
            dtg.dt,
            dte.classification.dt[, c(ref.column.name, "translational_regulation"), with = FALSE],
            by = ref.column.name,
            all.x = TRUE
        )

    } else {
        dtg.dt[, translational_regulation := NA]
    }

    if(nrow(annot.df) == 0){
        annot.df <- normalized.polysome.count.dt[, ref.column.name, with = FALSE]
    } else {"Use annot.df"}
    
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
