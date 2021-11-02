
library("magrittr")
library("data.table")

## Read differential translation analysis results
readDteResults <- function(dte.comparison.name, result.file.dir, analysis.level = "gene"){
    dte.comparison.filename <- dte.comparison.name %>%
        gsub("\\(", "", .) %>%
        gsub("\\)", "", .) %>%
        gsub("\\|", "-", .)
    
    dte.res.dt <- fread(file.path(
        result.file.dir,
        paste0(dte.comparison.filename, ".csv")
    ))
    dte.res.dt[, `:=`(
        analysis_level = analysis.level,
        comparison_name = dte.comparison.name
    )]
    
    return(dte.res.dt)
}

extractDteComparisonInfo <- function(dt){
    dt <- copy(dt)
    dt[, `:=`(
        cell = str_split_fixed(comparison_name, "_", n = 6)[, 1] %>%
            factor(., levels = c("RCC4", "786O")),
        VHL = str_split_fixed(comparison_name, "_", n = 6)[, 2] %>%
            factor(., levels = c("noVHL", "VHL")),
        EIF4E2 = str_split_fixed(comparison_name, "_", n = 6)[, 3] %>%
            factor(., levels = c("noEIF4E2", "EIF4E2")),
        clone = str_split_fixed(comparison_name, "_", n = 6)[, 4],
        treatment = str_split_fixed(comparison_name, "_", n = 6)[, 5],
        compared_elements = str_split_fixed(comparison_name, "__", n = 2)[, 2] %>%
            factor(., levels = c("noVHL_vs_VHL", "noEIF4E2_vs_EIF4E2", "Torin1_vs_DMSO")),
        plot_x_label = gsub("__", "\n", comparison_name) %>%
            factor(levels = unique(.)),
        experiment_purpose = case_when(
            comparison_name %in% c(
                                     "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL",
                                     "786O_xx_EIF4E2_yy_NA__noVHL_vs_VHL"
                                 ) ~ "VHL",
            comparison_name %in% c(
                                     "786O_xx_EIF4E2_1_NA__noVHL_vs_VHL",
                                     "786O_xx_noEIF4E2_1_NA__noVHL_vs_VHL",
                                     "786O_noVHL_xx_1_NA__noEIF4E2_vs_EIF4E2",
                                     "786O_VHL_xx_1_NA__noEIF4E2_vs_EIF4E2",
                                     "786O_(VHL|noVHL)_xx_1_NA__noEIF4E2_vs_EIF4E2"
                                 ) ~ "EIF4E2",
            comparison_name == "786O_noVHL_EIF4E2_1_xx__Torin1_vs_DMSO" ~ "mTOR"
        ) %>% factor(levels = c("VHL", "EIF4E2", "mTOR"))
    )]
    return(dt)
}

