

analyzeGoEnrichment <- function(ontology.type, input.genes, universe.genes){

    library("topGO")
    library("org.Hs.eg.db")
    
    input.gene.vec <- factor(as.integer(
        universe.genes %in% input.genes
    ))
    names(input.gene.vec) <- universe.genes
    GO.BP <- new(
        "topGOdata",
        ontology = ontology.type,
        allGenes = input.gene.vec,
        nodeSize = 10,
        annot = annFUN.org, 
        mapping = "org.Hs.eg.db",
        ID = "ENSEMBL"
    )
    GO.BP.wdight01.res <- runTest(
        GO.BP,
        algorithm = "weight01",
        statistic = "fisher"
    )
    GO.BP.parentchild.res <- runTest(
        GO.BP,
        algorithm = "parentchild",
        statistic = "fisher"
    )
    sig.GO.counts <- geneData(GO.BP.wdight01.res)["Significant"]
    go.analysis.out.dt <- GenTable(
        GO.BP,
        pvalue = GO.BP.wdight01.res,
        pc_p = GO.BP.parentchild.res,
        orderBy = "pvalue",
        topNodes = sig.GO.counts
    ) %>%
        data.table %>%
    cbind(
        ontology = ontology.type
    )
    return(go.analysis.out.dt)
}


analyzeGoEnrichments <- function(regulation.type, input.genes, universe.genes){

    library("topGO")
    library("org.Hs.eg.db")
        
    GO.all.res.dt <- suppressMessages(lapply(
        c("BP", "MF", "CC"),
        analyzeGoEnrichment,
        input.genes = input.genes,
        universe.genes = universe.genes
    )) %>%
        rbindlist %>%
    cbind(data.table(regulation_type = regulation.type))
    return(GO.all.res.dt)
}
