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
