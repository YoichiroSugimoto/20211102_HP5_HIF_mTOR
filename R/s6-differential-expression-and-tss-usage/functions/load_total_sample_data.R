## This script aims to load basic sample data files, and comparison parameters so that all the scripts in this directory becomes consistent

require("stringr")
require("data.table")

## Sample data preprocessing
sample.file <- normalizePath(list.files(
    "../../../",
    pattern = "^processed_sample_file.csv",
    recursive = TRUE,
    full.names = TRUE
))

print(paste0("Sample file used: ", sample.file))

sample.dt <- fread(sample.file)
sample.names <- sample.dt[, sample_name]
total.sample.dt <- sample.dt[experiment == "total"]

## Set clone so that really the same clone has the same number
total.sample.dt <- total.sample.dt[, `:=`(
    sample_name = gsub("^total_", "", sample_name) %>%
        {str_split_fixed(., "_", n = 6)[, 1:5]} %>% # remove technical replicate information
        data.table %>%
        {.[, do.call(paste, c(.SD, sep = "_"))]}
   ,
    clone =
        clone +
        ifelse(VHL == "VHL", 10, 0) +
        ifelse(HIF1B == "noHIF1B", 100, 0)
)][!duplicated(sample_name)]


## Prepare DESeq2 coldata
total.coldata.cols <- c("cell", "VHL", "HIF1B", "oxygen", "clone")
total.coldata.df <- data.frame(
    total.sample.dt[, total.coldata.cols, with = FALSE],
    row.names = total.sample.dt[, sample_name]
)
## Define the levels of factors
total.coldata.df[, "cell"] <- factor(
    total.coldata.df[, "cell"], levels = c("RCC4", "786O")
)
total.coldata.df[, "VHL"] <- factor(
    total.coldata.df[, "VHL"], levels = c("VHL", "noVHL")
)
total.coldata.df[, "oxygen"] <- factor(
    total.coldata.df[, "oxygen"], levels = c("N", "H")
)
total.coldata.df[, "HIF1B"] <- factor(
    total.coldata.df[, "HIF1B"], levels = c("HIF1B", "noHIF1B")
)
total.coldata.df[, "clone"] <- factor(
    total.coldata.df[, "clone"]
)


## List of comparisons
total.comparison.dt <- data.table(
    comparison = c(
        ## Deconvolute VHL dependent gene expression regulation
        "RCC4_xx_HIF1B_N__noVHL_vs_VHL",
        "786O_xx_HIF1B_N__noVHL_vs_VHL",
        "786O_xx_noHIF1B_N__noVHL_vs_VHL",
        ## Deconvolute hypoxia dependent gene expression regulation
        "RCC4_noVHL_HIF1B_xx__H_vs_N",
        "RCC4_VHL_HIF1B_xx__H_vs_N"
    ),
    exp_formula = c(
        ## Deconvolute VHL dependent gene expression regulation
        ~ VHL,
        ~ VHL,
        ~ VHL,
        ## Deconvolute hypoxia dependent gene expression regulation
        ~ clone + oxygen,
        ~ clone + oxygen
    )
)

print("The following R objects were exported: total.sample.dt, total.coldata.df, total.comparison.dt")
