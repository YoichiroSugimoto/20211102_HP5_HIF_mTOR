### This is a script to run an Rmd script from command line
### Usage: Rscript run_rmd.R rmd_filename.R

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(rmarkdown)
rmarkdown::render(args[1])

if(TRUE){
    full.rmd.path <- normalizePath(args[1])

    doc.dir <- file.path(
        gsub("/R/", "/doc/", dirname(full.rmd.path))
    )

    create.dir <- function(dir.name){
        if(dir.exists(dir.name) == FALSE) {
            dir.create(dir.name, recursive = TRUE)
        }
    }

    create.dir(doc.dir)

    ## Move knitr output files to doc dir
    non.rmd.files <- grep(
        list.files(
            dirname(full.rmd.path),
            pattern = gsub(".rmd", "", basename(full.rmd.path)),
            full.names = TRUE
        ),
        pattern = ".rmd$",
        inv = TRUE,
        value = TRUE
    )

    non.rmd.files <- basename(non.rmd.files[!file.info(non.rmd.files)$isdir])
    
    file.rename(
        from = file.path(dirname(full.rmd.path), non.rmd.files),
        to = file.path(doc.dir, non.rmd.files)
    )

    ## If figure directory exists
    non.rmd.dir <- file.path(gsub(".rmd", "_files", basename(full.rmd.path)))

    if(dir.exists(file.path(dirname(full.rmd.path), non.rmd.dir))){

        file.mv.cmd <- paste(
            "rsync", "--remove-source-files", "-r",
            file.path(dirname(full.rmd.path), non.rmd.dir),
            file.path(doc.dir)
        )

        system(file.mv.cmd)

        empty.file.rm.cmd <- paste(
            "find",
            file.path(dirname(full.rmd.path), non.rmd.dir),
            "-type", "d",
            "-empty", "-delete"
        )

        system(empty.file.rm.cmd)
        
    } else {"No output dir"}
    
} else {
    "knitr output not moved"
}
