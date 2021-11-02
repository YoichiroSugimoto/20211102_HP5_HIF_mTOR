tssWeightByPosition <- function(filtered.tss.with.quantile.dt, input.bam.files){
    
    filtered.tss.with.quantile.dt <- filtered.tss.with.quantile.dt[grepl("^ENSG", gene_id)]

    print("The following files will be used:")
    print(basename(input.bam.files))

    readAlltssBam <- function(bam.file){
        require("GenomicAlignments")
        bam <- readGAlignments(bam.file)
        bam.dt <- bam %>% as.data.frame %>% data.table
        bam.count.dt <- bam.dt[, list(count = .N), by = list(seqnames, strand, start, end)]
        return(bam.count.dt)
    }

    all.bam.long.dt <- mclapply(
        input.bam.files,
        readAlltssBam,
        mc.cores = processors
    ) %>%
        rbindlist

    bam.count.dt <- all.bam.long.dt[
      , list(count = sum(count)), by = list(seqnames, strand, start, end)
    ]

    ol.dt <- findOverlaps(
        makeGRangesFromDataFrame(bam.count.dt),
        makeGRangesFromDataFrame(filtered.tss.with.quantile.dt)
    ) %>% data.frame %>% data.table

    bam.count.tss.dt <- cbind(
        bam.count.dt[ol.dt[, queryHits]],
        filtered.tss.with.quantile.dt[
            ol.dt[, subjectHits],
            c("tss_name", "q0", "q25", "q50", "q75", "q100"),
            with = FALSE
        ]
    )

    bam.count.tss.dt[, tss_count := sum(count), by = tss_name]
    bam.count.tss.dt[, tss_idx := 1:.N, by = tss_name]
    bam.count.tss.dt[, tss_name_with_idx := paste0(tss_name, "_", tss_idx)]

    bam.count.tss.dt[
      , rel_tss_position := case_when(## Calculate the relative position from 0
            strand == "+" ~ start - q0 + 1,
            strand == "-" ~ q100 - start + 1
        )]

    bam.count.tss.dt <- bam.count.tss.dt[seqnames != "chrM"]
    
    return(bam.count.tss.dt)
}

extractTxMetaInfo <- function(tx.gtf.file, tss.count.per.pos.dt, ref.colname, tx.info.dir, tx.meta.data.file, processors){
    ## Export Tx sequences
    print("Export transcript sequences")
    tx.mean.len.list <- exportTxSeqs(
        tx.gtf.file = tx.gtf.file,
        tss.count.per.pos.dt = tss.count.per.pos.dt,
        ref.colname = ref.colname,
        tx.info.dir = tx.info.dir
    )

    tx.mean.len.dt <- tx.mean.len.list$tx.mean.len.dt
    all.tx.len.dt <- tx.mean.len.list$all.tx.len.dt

    print("Analyzing uORF")
    all.tx.seqs <- file.path(
        tx.info.dir, "sub-data",
        paste0("all_tx_seq", ".fa")
    ) %>% readDNAStringSet

    uorf.ref.count.dt <- uorfNumber(
        all.tx.seqs = all.tx.seqs,
        all.tx.len.dt = all.tx.len.dt,
        tss.count.per.pos.dt = tss.count.per.pos.dt,
        ref.colname = ref.colname,
        tx.mean.len.dt = tx.mean.len.dt
    )

    print("Analyzing TOP motif")
    first75mer.seqs.file <- file.path(
        tx.info.dir, "sub-data",
        "first_75_seq.fa"
    )

    top.len.per.tss.dt <- calculateTopLength(
        first75mer.seqs.file = first75mer.seqs.file,
        tss.count.per.pos.dt = tss.count.per.pos.dt,
        ref.colname = ref.colname
    )

    print("Analyzing RNA structure")
    all.rnalfold.dt <- calculateRNAStructure(
        tx.info.dir = tx.info.dir,
        tss.count.per.pos.dt = tss.count.per.pos.dt,
        ref.colname = ref.colname,
        processors = processors
    )

    print("Export all data")
    ## Export all outputs
    all.individual.dt.file <- file.path(
        tx.info.dir, "sub-data", "all-meta-data-dts.RData"
    )

    save(
        tss.count.per.pos.dt,
        tx.mean.len.dt,
        all.tx.len.dt,
        uorf.ref.count.dt,
        top.len.per.tss.dt,
        all.rnalfold.dt,
        file = all.individual.dt.file
    )

    all.tx.meta.dt <- Reduce(
        function(...) merge(..., all = TRUE, by = ref.colname),
        list(
            tx.mean.len.dt,
            uorf.ref.count.dt,
            top.len.per.tss.dt,
            all.rnalfold.dt
        )
    )

    uorf.cols <- c("uORF_frame0", "uORF_frame1", "uORF_frame2", "uORF_all")

    setnafill(
        all.tx.meta.dt,
        fill = 0,
        cols = uorf.cols
    )

    all.tx.meta.dt[, `:=`(
        tx_len = NULL,
        utr5_len = NULL,
        cds_MFE_G4p_local = NULL
    )]

    all.tx.meta.dt <- all.tx.meta.dt[
      , (uorf.cols) := lapply(.SD, function(x){if_else(is.na(min_utr5_len), NA_real_, x)}),
        .SDcols = uorf.cols
    ][!is.na(gene_id) & cds_len > 100 & !is.na(mean_utr5_len)]

    fwrite(all.tx.meta.dt, tx.meta.data.file)

    return()
}

## 1. Export Tx sequences
exportTxSeqs <- function(tx.gtf.file, tss.count.per.pos.dt, ref.colname, tx.info.dir){

    bs.genome.ver <- "BSgenome.Hsapiens.UCSC.hg38"
    library(bs.genome.ver, character.only = TRUE)

    ## Define directories
    tx.info.subdir <- file.path(tx.info.dir, "sub-data")
    create.dirs(c(tx.info.subdir))

    ## Tx sequence will be extracted
    promoter.assigned.txdb <- makeTxDbFromGFF(file = tx.gtf.file)

    tx.seqs <- extractTranscriptSeqs(
        BSgenome.Hsapiens.UCSC.hg38,
        promoter.assigned.txdb,
        use.names = TRUE
    )

    tx.len.dt <- data.table(transcriptLengths(
        promoter.assigned.txdb,
        with.utr5_len = TRUE,
        with.cds_len = TRUE,
        with.utr3_len = TRUE
    ))

    setnames(tx.len.dt, old = "tx_name", new = ref.colname)

    utr5.seqs <- narrow(tx.seqs, start = 1, end = tx.len.dt[, utr5_len])

    cds.seqs <- narrow(
        tx.seqs,
        start = tx.len.dt[, utr5_len] + 1L,
        end = -(tx.len.dt[, utr3_len] + 1L)
    ) 

    utr3.seqs <- narrow(
        tx.seqs,
        start = tx.len.dt[, tx_len] - tx.len.dt[, utr3_len] + 1
    ) 

    ## Export sequences
    tx.seqs.file <- file.path(
        tx.info.subdir,
        "transcript-sequence.fa"
    )
    export(tx.seqs, con = tx.seqs.file, format = "fasta")

    temp <- paste(
        "samtools", "faidx",
        tx.seqs.file,
        "-o", gsub(".fa$", ".fa.fai", tx.seqs.file)
    ) %>% system.cat
    
    cds.seqs.file <- file.path(
        tx.info.subdir,
        "cds-sequence.fa"
    )
    export(cds.seqs[width(cds.seqs) != 0], con = cds.seqs.file, format = "fasta")

    utr3.seqs.file <- file.path(
        tx.info.subdir,
        "utr3-sequence.fa"
    )
    export(utr3.seqs[width(utr3.seqs) != 0], con = utr3.seqs.file, format = "fasta")

    ## Export exact transcript sequence
    tx.seqs.dt <- tx.seqs %>% as.data.frame %>% data.table(keep.rownames = ref.colname)
    setnames(tx.seqs.dt, old = "x", new = "all_seq")

    all.tx.info.dt <- Reduce(
        function(...) merge(..., all = TRUE, by = ref.colname),
        list(
            tss.count.per.pos.dt,
            tx.len.dt,
            tx.seqs.dt
        )
    )

    all.tx.info.dt[, `:=`(
        ind_tx_len = (tx_len - rel_tss_position + 1) %>%
            {if_else(. < 100, NA_real_, .)}, # Ignore tx length shorter than  100 nts 
        ind_utr5_len = (utr5_len - rel_tss_position + 1) %>%
            {if_else(. < 0, NA_real_, .)}
    )]

    all.tx.info.dt[, `:=`(
        first_75_seq = str_sub(
            all_seq, start = rel_tss_position, end = rel_tss_position + 74
        ) %>% {if_else(ind_tx_len < 75, NA_character_, .)},
        non_first_75_utr5_seq = str_sub(
            all_seq,
            start = rel_tss_position + 75,
            end = rel_tss_position + ind_utr5_len - 1
        ),
        all_tx_seq = str_sub(
            all_seq,
            start = rel_tss_position
        )
    )]

    all.tx.info.dt[, `:=`(
        all_seq = NULL
    )]

    ## export sequences
    createNamedDNAStringSet <- function(seq.colname, dt, tx.info.subdir){
        dt <- dt[!is.na(get(seq.colname))]
        seqs <- DNAStringSet(dt[, get(seq.colname)])
        names(seqs) <- dt[, tss_name_with_idx]

        file.path(
            tx.info.subdir,
            paste0(seq.colname, ".fa")
        ) %>%
            {rtracklayer::export(
                              seqs[width(seqs) != 0],
                              con = .,
                              format = "fasta"
                          )}
        
        return(seqs)
    }

    temp <- lapply(
        c("all_tx_seq", "first_75_seq", "non_first_75_utr5_seq"),
        createNamedDNAStringSet,
        dt = all.tx.info.dt,
        tx.info.subdir = tx.info.subdir
    )

    temp <- file.path(
        tx.info.subdir,
        paste0("all_tx_seq", ".fa")
    ) %>% {paste(
               "samtools", "faidx",
               .,
               "-o", gsub(".fa$", ".fa.fai", .)
           )} %>% system.cat
    
    ind.len.dt <- all.tx.info.dt[
      , lapply(.SD, function(x){sum((count / tss_count) * x)}),
        .SDcols = c("ind_tx_len", "ind_utr5_len"),
        by = get(ref.colname)
    ]
    setnames(ind.len.dt, old = "get", new = ref.colname)

    min.utr5.len <- all.tx.info.dt[
      , list(min_utr5_len = min(ind_utr5_len)), by = get(ref.colname)
    ]
    setnames(min.utr5.len, old = "get", new = ref.colname)

    ## Calculate Kozak score
    library("ORFik")
    promoter.assigned.cds.gr <- cdsBy(
        promoter.assigned.txdb, by = "tx", use.names = TRUE
    )
    promoter.assigned.tx.gr <- exonsBy(
        promoter.assigned.txdb, by = "tx", use.names = TRUE
    )
    tx.kozak.scores <- kozakSequenceScore(
        promoter.assigned.cds.gr,
        promoter.assigned.tx.gr,
        BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
        species = "human"
    )
    tx.kozak.score.dt <- data.table(
        ref_col = names(promoter.assigned.cds.gr),
        tx_kozak_score = tx.kozak.scores 
    )
    setnames(tx.kozak.score.dt, old = "ref_col", new = ref.colname)

    tx.mean.len.dt <- Reduce(
        function(...) merge(..., all = TRUE, by = ref.colname),
        list(
            tx.len.dt[, c(
                ref.colname,
                "gene_id",
                "tx_len", "utr5_len", "cds_len", "utr3_len"
            ), with = FALSE],
            ind.len.dt,
            min.utr5.len,
            tx.kozak.score.dt
        )
    )

    setnames(
        tx.mean.len.dt,
        old = colnames(tx.mean.len.dt), new = gsub("ind_", "mean_", colnames(tx.mean.len.dt))
    )

    all.tx.len.dt <- all.tx.info.dt[, c(
        ref.colname, "tss_name_with_idx",
        "ind_tx_len", "ind_utr5_len", "cds_len", "utr3_len"
    ), with = FALSE]
    
    return(list(tx.mean.len.dt = tx.mean.len.dt, all.tx.len.dt = all.tx.len.dt))
}

## Calculate uORF number
uorfNumber <- function(all.tx.seqs, all.tx.len.dt, tss.count.per.pos.dt, ref.colname, tx.mean.len.dt){

    library("ORFik")

    tx.mRNA.seqs <- all.tx.seqs[
        names(all.tx.seqs) %in%
        all.tx.len.dt[
            get(ref.colname) %in% tx.mean.len.dt[!is.na(min_utr5_len), get(ref.colname)]
        ][cds_len > 0, tss_name_with_idx]
    ]

    all.orf.ranges <- ORFik::findORFs(
        tx.mRNA.seqs,
        startCodon = "ATG",
        longestORF = TRUE
    )

    names(all.orf.ranges) <- names(tx.mRNA.seqs)[as.integer(names(all.orf.ranges))]

    all.uorf.dt <- unlist(all.orf.ranges) %>%
        as.data.frame %>%
        data.table

    setnames(all.uorf.dt, old = "names", new = "tss_name_with_idx")

    all.uorf.dt <- merge(
        all.uorf.dt,
        all.tx.len.dt,
        by = "tss_name_with_idx",
        all.x = TRUE
    )

    all.uorf.dt[, `:=`(
        uORF_flag = start < ind_utr5_len + 1,
        frame = magrittr::mod(start - (ind_utr5_len + 1), 3),
        CDS_ol = if_else(end > ind_utr5_len, TRUE, FALSE),
        aa_len = (end - start - 2) / 3
    )]

    uorf.count.dt <- all.uorf.dt[uORF_flag == TRUE][
      , list(uORF_count = .N), by = list(tss_name_with_idx, frame)
    ]

    d.uorf.count.dt <- dcast(
        uorf.count.dt,
        tss_name_with_idx ~ frame,
        value.var = "uORF_count"
    )

    setnafill(d.uorf.count.dt, fill = 0, cols = c("1", "2", "0"))

    setnames(
        d.uorf.count.dt,
        old = c("0", "1", "2"),
        new = c("uORF_frame0", "uORF_frame1", "uORF_frame2")
    )

    d.uorf.count.dt[, `:=`(
        uORF_all = uORF_frame0 + uORF_frame1 + uORF_frame2
    )]

    d.uorf.count.dt <- merge(
        d.uorf.count.dt,
        tss.count.per.pos.dt[, c(
            "tss_name_with_idx", ref.colname,
            "count", "tss_count"
        ), with = FALSE],
        by = "tss_name_with_idx"
    )

    uorf.ref.count.dt <- d.uorf.count.dt[
      , lapply(.SD, function(x){sum((count / tss_count) * x)}),
        .SDcols = grep("^uORF", colnames(d.uorf.count.dt), value = TRUE),
        by = get(ref.colname)
    ]

    setnames(uorf.ref.count.dt, old = "get", new = ref.colname)

    return(uorf.ref.count.dt)
}

calculateTopLength <- function(first75mer.seqs.file, tss.count.per.pos.dt, ref.colname){

    first75mer.seqs <- readDNAStringSet(first75mer.seqs.file)

    first75me.seqs.dt <- first75mer.seqs %>% as.data.frame %>%
        data.table(keep.rownames = "tss_name_with_idx")
    setnames(first75me.seqs.dt, old = "x", new = "first_75_mer")
    
    ## first75me.seqs.dt should have first_75_mer column
    first75me.seqs.dt[, `:=`(
        first_C_position = str_locate(first_75_mer, "C")[, "end"],
        first_T_position = str_locate(first_75_mer, "T")[, "end"],
        first_p_position = str_locate(first_75_mer, "(C|T)")[, "end"],
        TOP_length =
            str_locate(first_75_mer, "C(T|C)*") %>% {.[, "end"] - .[, "start"] + 1} %>%
            {pmin(., 10)},
        tTOP_length =
            str_locate(first_75_mer, "T(T|C)*") %>% {.[, "end"] - .[, "start"] + 1} %>%
            {pmin(., 10)},
        pTOP_length =
            str_locate(first_75_mer, "(T|C)(T|C)*") %>% {.[, "end"] - .[, "start"] + 1} %>%
            {pmin(., 10)},
        is_first_base_G = grepl("^G", first_75_mer),
        PRTE_score =
            case_when(
                !grepl("^(A|T|C|G){5}T", first_75_mer) ~ as.integer(0),
                TRUE ~ str_count(substr(first_75_mer, start = 1, stop = 9), "(T|C)")
            )
    )]

    first75me.seqs.dt[, `:=`(
        p1_TOP = if_else(first_C_position == 1, TOP_length, 0, missing = 0),
        p2_TOP = if_else(first_C_position == 2, TOP_length, 0, missing = 0),
        p3_TOP = if_else(first_C_position == 3, TOP_length, 0, missing = 0),
        p1_tTOP = if_else(first_T_position == 1, tTOP_length, 0, missing = 0),
        p2_tTOP = if_else(first_T_position == 2, tTOP_length, 0, missing = 0),
        p3_tTOP = if_else(first_T_position == 3, tTOP_length, 0, missing = 0),
        p1_pTOP = if_else(first_p_position == 1, pTOP_length, 0, missing = 0),
        p2_pTOP = if_else(first_p_position == 2, pTOP_length, 0, missing = 0),
        p3_pTOP = if_else(first_p_position == 3, pTOP_length, 0, missing = 0)
    )]

    first75me.seqs.dt[, `:=`(
        first_C_position = NULL,
        first_T_position = NULL,
        first_p_position = NULL,
        TOP_length = NULL,
        tTOP_length = NULL,
        pTOP_length = NULL
    )]

    top.len.dt <- merge(
        tss.count.per.pos.dt,
        first75me.seqs.dt,
        by = "tss_name_with_idx"
    )

    top.len.per.tss.dt <- top.len.dt[
      , lapply(.SD, function(x){sum((count / tss_count) * x)}),
        .SDcols = grep("TOP$|PRTE_score$", colnames(top.len.dt), value = TRUE),
        by = get(ref.colname)
    ]
    setnames(top.len.per.tss.dt, old = "get", new = ref.colname)

    grep("TOP$|PRTE_score$", colnames(top.len.per.tss.dt), value = TRUE) %>%
        {setnames(top.len.per.tss.dt, old = ., new = paste0("tss_", .))}

    return(top.len.per.tss.dt)
}

calculateRNAStructure <- function(tx.info.dir, tss.count.per.pos.dt, ref.colname, processors){
    
    runRNALfold <- function(input.fa.file, mfe.name, additional.args = "", ref.colname = "tss_name"){

        rnalfold.out.file <- file.path(
            dirname(input.fa.file),
            paste0(mfe.name, "_RNALfold.txt")
        )

        rnalfold.cmd <- paste(
            "RNALfold",
            additional.args,
            ## Help message was not correct!!
            "<", basename(input.fa.file),
            ">", rnalfold.out.file
        )

        system(rnalfold.cmd)

        rnalfold.out <- readLines(rnalfold.out.file)

        splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

        split.lfold.outs <- splitAt(rnalfold.out, grep(">ENS", rnalfold.out))

        extractMFE <- function(rnalfold.out){
            local.mfe <- str_split_fixed(rnalfold.out[2:(length(rnalfold.out) - 2)], " \\(", n = 2)[, 2]

            local.mfe.dt <- data.table(
                tss_name = gsub(">", "", rnalfold.out[1]),
                structure_start = str_split_fixed(local.mfe, "\\)", n = 2)[, 2] %>%
                    gsub(" ", "", .) %>%
                    as.numeric(.),
                mfe____local = str_split_fixed(local.mfe, "\\)", n = 2)[, 1] %>%
                    gsub(" ", "", .) %>%
                    as.numeric(.),
                mfe___ = as.numeric(gsub("( \\(|\\))", "", rnalfold.out[length(rnalfold.out)]))
            )

            return(local.mfe.dt)
        }

        lfold.dt <- lapply(
            split.lfold.outs,
            extractMFE
        ) %>%
            rbindlist

        setnames(lfold.dt, old = "tss_name", new = ref.colname)
        
        setnames(
            lfold.dt,
            old = colnames(lfold.dt),
            new = gsub("mfe___", mfe.name, colnames(lfold.dt))
        )

        mfe.out.file <- file.path(
            dirname(input.fa.file),
            paste0(mfe.name, "_long.csv")
        )
        
        fwrite(lfold.dt, mfe.out.file)

        lfold.summary.dt <- lfold.dt
        
        lfold.summary.dt[, min_local := min(get(paste0(mfe.name, "_local"))), by = get(ref.colname)]
        
        lfold.summary.dt <- lfold.summary.dt[
            get(paste0(mfe.name, "_local")) == min_local
        ][!duplicated(get(ref.colname))]

        lfold.summary.dt[, `:=`(
            structure_start = NULL,
            min_local = NULL
        )]
        
        mfe.out.summary.file <- file.path(
            dirname(input.fa.file),
            paste0(mfe.name, ".csv")
        )
        
        fwrite(lfold.summary.dt, mfe.out.summary.file)
        
        return()
    }


    tx.segment.seq.files <- c("first_75_seq.fa", "non_first_75_utr5_seq.fa", "cds-sequence.fa") %>%
        {file.path(tx.info.dir, "sub-data", .)} %>%
        lapply(normalizePath)

    tx.segments <- c("first75mer", "nonFirst75mer", "cds")


    ## This is to run RNALfold
    curwd <- getwd()
    setwd(file.path(tx.info.dir, "sub-data"))

    temp <- mcmapply(
        runRNALfold,
        tx.segment.seq.files,
        paste0(tx.segments, "_MFE_G4p"),
        additional.args = "--gquad",
        ref.colname = c("tss_name_with_idx", "tss_name_with_idx", ref.colname),
        mc.cores = processors
    )

    setwd(curwd)

    cds.rnalfold.dt <- file.path(
        tx.info.dir, "sub-data", "cds_MFE_G4p.csv"
    ) %>%
        fread

    calculateTssWeightedMFE <- function(structure.name, tx.info.dir, tss.count.per.pos.dt, ref.colname){

        top75.rnalfold.long.dt <- file.path(
            tx.info.dir, "sub-data",
            paste0(
                structure.name, "_MFE_G4p_long.csv"
            )
        ) %>%
            fread
        
        top75.rnalfold.long.dt[, min_local := min(
                                     get(paste0(structure.name, "_MFE_G4p_local"))
                                 ), by = tss_name_with_idx]


        top75.rna.dt <- top75.rnalfold.long.dt[
            get(paste0(structure.name, "_MFE_G4p_local")) == min_local
        ][!duplicated(tss_name_with_idx)]

        top75.rna.structure.sum.dt <- merge(
            tss.count.per.pos.dt[
              , c("tss_name_with_idx", ref.colname, "count", "tss_count"), with = FALSE
            ],
            top75.rna.dt,
            by = "tss_name_with_idx"
        )
        
        final.top75.rna.structure.sum.dt <- top75.rna.structure.sum.dt[
          , lapply(.SD, function(x){sum((count / tss_count) * x)}),
            .SDcols = paste0(structure.name, "_MFE_G4p"),
            by = get(ref.colname)
        ]

        setnames(final.top75.rna.structure.sum.dt, old = "get", new = ref.colname)
        
        return(final.top75.rna.structure.sum.dt)
    }

    utr5.rnalfold.dt <- lapply(
        c("first75mer", "nonFirst75mer"),
        calculateTssWeightedMFE,
        tx.info.dir = tx.info.dir,
        tss.count.per.pos.dt = tss.count.per.pos.dt,
        ref.colname = ref.colname
    ) %>%
        {Reduce(
             function(...) merge(..., all = TRUE, by = ref.colname), .
         )}

    all.rnalfold.dt <- merge(
        utr5.rnalfold.dt,
        cds.rnalfold.dt,
        by = ref.colname,
        all = TRUE
    )

    return(all.rnalfold.dt)
}

