## This function load the analysis results by s6-2 this includes:
## 1. TSS level differential expression
## 2. TSS level differential usage ratio
## 3. Differential TSS analysis by combining 1 and 2

library("data.table")
library("stringr")

comparison.info.file <- normalizePath(list.files(
    file.path("../"),
    pattern = "load_total_sample_data.R$",
    recursive = TRUE,
    full.names = TRUE
))

comparison.info.file <- grep(
    "s6-differential-expression-and-tss-usage/functions",
    comparison.info.file, value = TRUE
)

source(comparison.info.file)
print("Comparison information was loaded")

all.dirs <- normalizePath(list.dirs(
    file.path("../../../")
))

results.dir <- grep("/results$", all.dirs, value = TRUE)

print(results.dir)

if(length(results.dir) != 1){
    stop("Unique results dir not found")
} else {"Found"}

s6.dir <- file.path(results.dir, "s6-differential-regulation-analysis")
s6.2.dir <- file.path(s6.dir, "s6-2-differentially-TSS")
s6.2.1.tss.ratio.dir <- file.path(s6.2.dir, "TSS_ratio_change")
s6.2.2.tss.de.dir <- file.path(s6.2.dir, "TSS_differential_expression")
s6.2.3.diff.tss.dir <- file.path(s6.2.dir, "Differential_TSS_usage")

## TSS information
s4.tss.dir <- file.path(results.dir, "s4-tss-definition-and-tx-assignment")
s4.1.tss.def.dir <- file.path(s4.tss.dir, "s4-1-tss-definition")
s4.1.6.filtered.tss.dir <- file.path(s4.1.tss.def.dir, "s4-1-6-filtered-tss")
filtered.tss.with.quantile.file <- file.path(
    s4.1.6.filtered.tss.dir,
    "filtered-tss-with-quantile.csv"
)

filtered.tss.with.quantile.dt <- fread(filtered.tss.with.quantile.file)

## Import analysis results
extractComparisonInfo <- function(input.dt){
    sl.dt <- copy(input.dt)
    
    sl.dt[, `:=`(
        comparison_base = str_split_fixed(comparison_name, "__", n = 2)[, 1] %>%
            gsub("log2FC_", "", .) %>%
            gsub("_xx", "", .),
        compared_elements = factor(
            str_split_fixed(comparison_name, "__", n = 2)[, 2],
            levels = c("noVHL_vs_VHL", "H_vs_N")),
        plot_label = paste0(gsub("__", "\n(", comparison_name), ")")
    )]

    sl.dt[, `:=`(
        comparison_base = factor(comparison_base, levels = unique(comparison_base)),
        plot_label = factor(plot_label, levels = unique(plot_label))
    )]

    return(sl.dt)
}

## 1. TSS level differential expression
readTssDeRes <- function(comparison.name, s6.2.dir){
    sl.dt <- fread(file.path(s6.2.dir, paste0(comparison.name, "_DE.csv")))
    sl.dt[, comparison_name := comparison.name]
    return(sl.dt)
}

tss.de.res.dt <- readTssDeRes(
    "RCC4_xx_HIF1B_N__noVHL_vs_VHL",
    s6.2.2.tss.de.dir
)

tss.de.res.dt <- extractComparisonInfo(tss.de.res.dt)

tss.de.res.dt <- merge(
    filtered.tss.with.quantile.dt[, .(tss_name, annot)],
    tss.de.res.dt,
    by = "tss_name"
)


## 2. TSS level differential usage ratio
readRatioRes <- function(comparison.name, s6.2.1.tss.ratio.dir){
    sl.dt <- fread(file.path(s6.2.1.tss.ratio.dir, paste0(comparison.name, ".csv")))
    sl.dt[, comparison_name := comparison.name]
    return(sl.dt)
}

tss.ratio.res.dts <- lapply(
    total.comparison.dt[, comparison],
    readRatioRes,
    s6.2.1.tss.ratio.dir = s6.2.1.tss.ratio.dir
)

tss.ratio.res.dt <- rbindlist(tss.ratio.res.dts, fill = TRUE)
tss.ratio.res.dt <- extractComparisonInfo(tss.ratio.res.dt)

## 3. Differential TSS analysis results
readDiffTssRes <- function(comparison.name, s6.3.dir){
    csv.file <- file.path(s6.3.dir, paste0(comparison.name, "-diff-TSS.csv"))
    if(file.exists(csv.file)){
        sl.dt <- fread(file.path(s6.3.dir, paste0(comparison.name, "-diff-TSS.csv")))
        sl.dt[, comparison_name := comparison.name]
    } else {
        sl.dt <- data.table()
    }
    return(sl.dt)
}

diff.tss.res.dt <- readDiffTssRes(
    "RCC4_xx_HIF1B_N__noVHL_vs_VHL",
    s6.3.dir = s6.2.3.diff.tss.dir 
)
diff.tss.res.dt <- extractComparisonInfo(diff.tss.res.dt)

print("The following objects were loaded: tss.de.res.dt, tss.ratio.res.dt, diff.tss.res.dt")
