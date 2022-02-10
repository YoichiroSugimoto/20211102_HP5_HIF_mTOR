trslWideToLong <- function(all.de.dte.res.dt, sl.comp.names, intervention.names){
    require("data.table")
    require("stringr")
    
    comp2int <- setNames(
        intervention.names,
        nm = sl.comp.names
    )

    m.all.de.dte.res.dt <- melt(
        all.de.dte.res.dt,
        id.vars = "gene_id",
        measure.vars = paste0(
            "MRL_log2fc_",
            sl.comp.names
        ),
        value.name = "MRL_log2FC"
    ) %>%
        {.[, `:=`(
             Intervention = comp2int[variable] %>%
                 factor(levels = rev(intervention.names))
         )]}
    
    trsl.cols <- paste0("translational_regulation_", sl.comp.names)
    trsl.annot.dt <- copy(all.de.dte.res.dt)

    trsl.annot.dt[, (trsl.cols) := lapply(
                        .SD,
                        function(x){
                            case_when(
                                x %in% c("Up", "Down") ~ as.character(x),
                                TRUE ~ "Not significant"
                            )
                        }
                    ), .SDcols = trsl.cols]

    m.trsl.annot.dt <- melt(
        trsl.annot.dt,
        id.vars = "gene_id",
        measure.vars = trsl.cols,
        value.name = "translation_regulation"
    ) %>%
        {.[, `:=`(
             Intervention = comp2int[variable] %>%
                 factor(levels = rev(intervention.names)),
             translation_regulation = factor(
                 translation_regulation,
                 levels = c("Not significant", "Down", "Up")
             )
         )]}

    m.all.de.dte.res.dt <- merge(
        m.all.de.dte.res.dt,
        m.trsl.annot.dt[, .(gene_id, Intervention, translation_regulation)],
        by = c("gene_id", "Intervention")
    ) %>% {.[, Intervention := factor(Intervention, levels = intervention.names)]}

    return(m.all.de.dte.res.dt)
}

plotTrslDistByIntervention <- function(m.all.de.dte.res.dt, show.quantile = TRUE){
    require("ggplot2")
    g1 <- ggplot(
            data = m.all.de.dte.res.dt,
        ) +
        geom_vline(xintercept = 0, color = "black") +
        facet_grid(Intervention ~ .) +
        coord_cartesian(xlim = c(-1.0, 0.5)) +
        theme(
            aspect.ratio = 0.5
        ) +
        ylab("Proportion") +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))

    if(show.quantile){
        q.dt <- m.all.de.dte.res.dt[, list(
            q25 = quantile(MRL_log2FC, probs = 0.25, na.rm = TRUE),
            q50 = quantile(MRL_log2FC, probs = 0.5, na.rm = TRUE),
            q75 = quantile(MRL_log2FC, probs = 0.75, na.rm = TRUE),
            ymin = -0,
            ymax = Inf
        ), by = Intervention]
        g1 <- g1 +
            geom_rect(
                data = q.dt,
                aes(
                    xmin = q25, xmax = q75,
                    ymin = ymin, ymax = ymax
                ),
                col = NA,
                fill = "red",
                alpha = 0.2,
                inherit.aes = FALSE
            )
    } else {"OK"}
    g1 <- g1 +
        geom_histogram(
            aes(
                x = MRL_log2FC,
                y = stat(count / sum(count))
            ),
            binwidth = 0.025
        )
    print(g1)
    
    return()
}


proportionPerFraction <- function(count.per.fraction.dt, ref.col = "gene_id", data.col.grep = "^RCC4"){

    base.cols <- c("gene_id", "gene_name", "biotype")

    m.count.per.fraction.dt <- melt(
        count.per.fraction.dt,
        id.vars = unique(c(ref.col, base.cols)),
        measure.vars = grep(data.col.grep, colnames(count.per.fraction.dt), value = TRUE),
        value.name = "normalized_count",
        variable.name = "sample_name"
    )

    m.count.per.fraction.dt[, `:=`(
        sample_group = sub("(.*?_.*?_.*?)_.*", "\\1", sample_name),
        cell = str_split_fixed(sample_name, "_", n = 8)[, 1],
        VHL = str_split_fixed(sample_name, "_", n = 8)[, 2],
        EIF4E2 = str_split_fixed(sample_name, "_", n = 8)[, 3],
        clone = str_split_fixed(sample_name, "_", n = 8)[, 5],
        treatment = str_split_fixed(sample_name, "_", n = 8)[, 6],
        fraction = str_split_fixed(sample_name, "_", n = 8)[, 7]
    )]

    m.count.per.fraction.dt[, `:=`(
        sample_group = paste(cell, VHL, EIF4E2, treatment, sep = "_")
    )]

    ## Across fractions
    m.count.per.fraction.dt[, `:=`(
        sum_across_fraction = sum(normalized_count, na.rm = TRUE)
    ), by = list(get(ref.col), gene_id, gene_name, sample_group, clone)]

    m.count.per.fraction.dt[, `:=`(
        ratio_across_fraction = normalized_count / sum_across_fraction
    )]

    m.count.per.fraction.dt[, `:=`(
        sd_norm_count = sd(normalized_count),
        mean_norm_count = mean(normalized_count),
        mean_ratio = mean(ratio_across_fraction),
        sd_ratio = sd(ratio_across_fraction),
        lower_ratio_range = mean(ratio_across_fraction) - sd(ratio_across_fraction),
        upper_ratio_range = mean(ratio_across_fraction) + sd(ratio_across_fraction)
    ), by = list(get(ref.col), gene_id, gene_name, sample_group, fraction)]

    ## By fractions (only useful for TSS data)
    m.count.per.fraction.dt[, `:=`(
        sum_by_fraction = sum(normalized_count, na.rm = TRUE)
    ), by = list(gene_id, gene_name, sample_group, clone, fraction)]

    m.count.per.fraction.dt[, `:=`(
        ratio_by_fraction = normalized_count / sum_by_fraction
    )]

    m.count.per.fraction.dt[, `:=`(
        mean_ratio_by_fraction = mean(ratio_by_fraction),
        sd_ratio_by_fraction = sd(ratio_by_fraction),
        lower_ratio_range_by_fraction = mean(ratio_by_fraction) - sd(ratio_by_fraction),
        upper_ratio_range_by_fraction = mean(ratio_by_fraction) + sd(ratio_by_fraction)
    ), by = list(get(ref.col), gene_id, gene_name, sample_group, fraction)]

    m.count.per.fraction.dt <- m.count.per.fraction.dt[clone == 3]

    return(m.count.per.fraction.dt)
}