## collapseTechnicalReplicates <- function(tss.count.dt, ref.name = "tss_name", technical.replicate.index.pos = 7){
##     ## This function assume that the technical replicates are indicated by the last part of the sample_name (e.g. ..._1, or ..._2)
    
##     m.tss.count.dt <- melt(
##         tss.count.dt,
##         id.vars = ref.name,
##         value.name = "count",
##         variable.name = "sample_name_with_tr"
##     ) %>%
##         {.[, sample_name := str_split_fixed(
##                  sample_name_with_tr,
##                  "_",
##                  n = technical.replicate.index.pos
##              )[, 1:(technical.replicate.index.pos - 1)] %>%
##                  data.table %>%
##              {.[, do.call(paste, c(.SD, sep = "_"))]}
##            ]}
    
##     cm.tss.count.dt <- m.tss.count.dt[
##       , list(collapsed_count = sum(count)), by = list(sample_name, get(ref.name))
##     ]
##     setnames(cm.tss.count.dt, old = "get", new = ref.name)
##     collapsed.tss.count.dt <- dcast(
##         cm.tss.count.dt,
##         paste0(ref.name, " ~ sample_name"),
##         value.var = "collapsed_count"
##     )
##     return(collapsed.tss.count.dt)
## }
