## This script aims to load sample data for polysome profiling experiments, and comparison parameters so that all the scripts in this directory becomes consistent
require("data.table")

## Sample data preprocessing
sample.file <- normalizePath(list.files(
    "../../../",
    pattern = "^processed_sample_file.csv$",
    recursive = TRUE,
    full.names = TRUE
))

print(paste0("Sample file used: ", sample.file))

if(length(sample.file) > 1){
    stop("Multiple sample files found...")
} else {"OK"}

sample.dt <- fread(sample.file)

poly.sample.dt <- sample.dt[experiment == "polysome"]
poly.sample.dt[is.na(poly.sample.dt)] <- "NA"


## Set clone so that really the same clone has the same number
poly.sample.dt[, clone :=
                     clone +
                     ifelse(VHL == "VHL", 10, 0)
               ]


poly.sample.dt[
   , sample_name := gsub("polysome_", "", sample_name)
    ]

## Prepare DESeq2 poly.coldata
poly.coldata.cols <- c(
    "cell", "VHL", "EIF4E2", "gRNA_id", "clone", "treatment", "fraction"
)
poly.coldata.df <- data.frame(
    poly.sample.dt[, poly.coldata.cols, with = FALSE],
    row.names = poly.sample.dt[, sample_name]
)

## Define the levels of factors
poly.coldata.df[, "cell"] <- factor(
    poly.coldata.df[, "cell"], levels = c("RCC4", "786O")
)
poly.coldata.df[, "VHL"] <- factor(
    poly.coldata.df[, "VHL"], levels = c("VHL", "noVHL")
)
poly.coldata.df[, "EIF4E2"] <- factor(
    poly.coldata.df[, "EIF4E2"], levels = c("EIF4E2", "noEIF4E2")
)
poly.coldata.df[, "clone"] <- factor(
    poly.coldata.df[, "clone"]
)
poly.coldata.df[, "treatment"] <- factor(
    poly.coldata.df[, "treatment"], levels = c("NA", "DMSO", "Torin1")
)

## List of comparisons
translation.comparison.dt <- data.table(
    comparison = c(
        ## Comparison includes all the information except guide RNA information where the wildcard will be used all the cases
        ## Deconvolute VHL/mTOR dependent translational regulation
        "RCC4_xx_EIF4E2_yy_NA__noVHL_vs_VHL",
        "RCC4_noVHL_EIF4E2_yy_xx__Torin1_vs_NA",
        "RCC4_VHL_EIF4E2_yy_xx__Torin1_vs_NA",
        "786O_xx_EIF4E2_yy_NA__noVHL_vs_VHL",
        ## Deconvolute EIF4E2 dependent translational regulation
        "786O_noVHL_xx_yy_NA__noEIF4E2_vs_EIF4E2",
        "786O_VHL_xx_yy_NA__noEIF4E2_vs_EIF4E2",
        "786O_(VHL|noVHL)_xx_yy_NA__noEIF4E2_vs_EIF4E2"
    ),
    exp_formula = c(
        ## Deconvolute VHL dependent translational regulation
        ~ fraction + VHL + fraction:VHL,
        ~ fraction + treatment + fraction:treatment,
        ~ fraction + treatment + fraction:treatment,
        ~ fraction + VHL + fraction:VHL,
        ## Deconvolute EIF4E2 dependent translational regulation
        ~ fraction + EIF4E2 + fraction:EIF4E2,
        ~ fraction + EIF4E2 + fraction:EIF4E2,
        ~ fraction + VHL + EIF4E2 + fraction:VHL + fraction:EIF4E2
    )
)

rm(
    poly.coldata.cols, sample.file, sample.dt
)

print("The following objects are exported: poly.coldata.df, poly.sample.dt, translation.comparison.dt")
print("In translation.comparison.dt, xx specifies the factor compared where the comparison is specified after __, while yy is a wildcard. From left, each factor specifies cell, VHL, EIF4E2, clone, and treatment")


