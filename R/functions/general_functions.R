## Some useful functions

create.dir <- function(dir.name){
    if(dir.exists(dir.name) == FALSE) {
        dir.create(dir.name)
    }
}

create.dirs <- function(dirs){
    for(dir.name in dirs){
        create.dir(dir.name)
    }
}

system.cat <- function(cmd){
    stdout.text <- system(cmd, intern = TRUE)
    return(stdout.text)
}

dropColumnDt <- function(dt, drop.vec){
    dt[ , !names(dt) %in% c(drop.vec), with = FALSE]
}

dropColumnDf <- function(df, drop.vec){
    df[ , !names(df) %in% c(drop.vec)]
}

exportSourceData <- function(
                             dt, original.colnames, export.colnames,
                             export.file.name
                             ){
    export.dt <- copy(dt)
    
    setnames(
        export.dt,
        old = original.colnames,
        new = export.colnames
    )
    fwrite(
        export.dt[, export.colnames, with = FALSE],
        file = file.path(source.data.by.panel.dir, export.file.name)
    )
    return()
}

## merge.dt.list <- function(dt.list){
##     ## https://gist.github.com/reinholdsson/67008ee3e671ff23b568
##     merged <- Reduce(function(...) merge(..., all = T), dt.list)
##     return(merged)
## }

