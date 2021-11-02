assignTxToTss <- function(tx.map.out.dir, filtered.tss.with.quantile.dt, all.tx.gtf, all.primary.tx.dt, sl.sample.names, star.aligned.bam.dir, software.path.dt, cell.name = "RCC4", function.dir){

    indiv.bams.dir <- file.path(tx.map.out.dir, "individual-bams")
    stringtie.in.dir <- file.path(tx.map.out.dir, "stringtie-in")
    stringtie.out.dir <- file.path(tx.map.out.dir, "stringtie-out")

    create.dirs(c(
        indiv.bams.dir,
        stringtie.in.dir,
        stringtie.out.dir
    ))

    tx.per.tss.gtf <- file.path(
        tx.map.out.dir,
        paste0("transcripts-per-TSS-for-", cell.name, ".gtf")
    )


    sample.groups <- unique(paste(
        str_split_fixed(sl.sample.names, "_", n = 6)[, 2],
        str_split_fixed(sl.sample.names, "_", n = 6)[, 3],
        str_split_fixed(sl.sample.names, "_", n = 6)[, 4],
        str_split_fixed(sl.sample.names, "_", n = 6)[, 5],
        sep = "_"
    ))

    print("The following sample groups are examined:")
    print(sample.groups)

    print("Bam files are split by TSS")
    promoters.dt <- splitBamByTss(
        sl.sample.names = sl.sample.names,
        filtered.tss.with.quantile.dt = filtered.tss.with.quantile.dt,
        tx.map.out.dir = tx.map.out.dir,
        star.aligned.bam.dir = star.aligned.bam.dir,
        indiv.bams.dir = indiv.bams.dir,
        function.dir = function.dir
    )

    for(sl.tss.index in 1:max(promoters.dt[, tss_index])){
        file.path(
            stringtie.in.dir,
            paste0("promoters_idx_", sl.tss.index, ".txt")
        ) %>%
            {fwrite(
                 promoters.dt[tss_index == sl.tss.index, list(
                                                             chr,
                                                             q50,
                                                             strand,
                                                             tss = "TSS"
                                                         )],
                 .,
                 col.names = FALSE, sep = "\t"
             )}
    }

    ## add XS tag to bam file to run stringtie
    temp <- lapply(
        1:max(promoters.dt[, tss_index]),
        addXStag2Bam,
        indiv.bams.dir = indiv.bams.dir,
        stringtie.in.dir = stringtie.in.dir,
        tagXSstrandedData.path =
            software.path.dt[software_name == "tagXSstrandedData", path]
    )

    print("Run StringTie")
    temp <- lapply(
        1:max(promoters.dt[, tss_index]),
        runStringTie,
        stringtie.in.dir = stringtie.in.dir,
        stringtie.out.dir = stringtie.out.dir,
        all.tx.gtf = all.tx.gtf
    )

    print("Assign Tx to TSS")
    ref.tx.assign.list <- refStringTieAnnotMatching(
        promoters.dt = promoters.dt,
        all.tx.gtf = all.tx.gtf,
        all.primary.tx.dt = all.primary.tx.dt,
        stringtie.out.dir = stringtie.out.dir
    )

    tx.assigned.tss.dt <- ref.tx.assign.list$tx.assigned.tss.dt
    ref.dt <- ref.tx.assign.list$ref.dt
    query.exon.dt = ref.tx.assign.list$query.exon.dt

    print("Export GTF for assigned Tx")
    tss.redef.ref.dt <- createGtfForAssignedTx(
        tx.assigned.tss.dt = tx.assigned.tss.dt,
        ref.dt = ref.dt,
        tx.per.tss.gtf = tx.per.tss.gtf,
        query.exon.dt = query.exon.dt
    )

    print("Exporting data")
    assigned.tx.info.dt <- tss.redef.ref.dt[!duplicated(tss_name)][, .(
                                               tss_name, seqnames, start, end, width, strand,
                                               original_transcript_id, CDS_canonical
                                           )]

    setnames(
        assigned.tx.info.dt,
        old = "original_transcript_id",
        new = "transcript_id"
    )


    assigned.tx.info.dt <- merge(
        assigned.tx.info.dt,
        all.primary.tx.dt,
        by = "transcript_id"
    )

    tss.assigned.tx.info.file <- file.path(
        tx.map.out.dir,
        paste0("TSS-assigned-tx-info-for-", cell.name, ".csv")
    )

    fwrite(
        assigned.tx.info.dt,
        tss.assigned.tx.info.file
    )

    return()
}


splitBamByTss <- function(sl.sample.names, filtered.tss.with.quantile.dt, tx.map.out.dir, star.aligned.bam.dir, indiv.bams.dir, function.dir){
    
    promoters.dt <- filtered.tss.with.quantile.dt[
        annot != "intergenic" ## Process intronic TSS in a different manner
    ]

    promoters.dt <- promoters.dt[
        order(gene_id, chr, ifelse(strand == "+", start, -1 * start))
    ]
    promoters.dt[, `:=`(
        tss_index = str_split_fixed(tss_name, "_", n = 2)[, 2] %>% as.integer,
        idx = 1:nrow(promoters.dt)
    )]

    promoters.file <- file.path(
        tx.map.out.dir,
        "promoter-summary.csv"
    )

    fwrite(promoters.dt, promoters.file)
    
    all.bam.inputs <- list.files(
        star.aligned.bam.dir,
        pattern = paste0(
            "(",
            paste(paste0(sl.sample.names, ".bam"), collapse = "|"),
            ")"
        ),
        full.names = TRUE
    )

    runSplitBamByTssIdx.cmd <- function(bam.file, promoters.file, indiv.bams.dir, function.dir){

        runSplitBamByTssIdx.rscript <- file.path(function.dir, "splitBamByTssIdx.R")

        cmd <- paste(
            "Rscript",
            runSplitBamByTssIdx.rscript,
            "-b", bam.file,
            "-p", promoters.file,
            "-o", indiv.bams.dir
        )

        system.cat(cmd)
        return()
    }

    temp <- mclapply(
        all.bam.inputs,
        runSplitBamByTssIdx.cmd,
        promoters.file = promoters.file,
        indiv.bams.dir = indiv.bams.dir,
        function.dir = function.dir,
        mc.cores = processors
    )

    return(promoters.dt)
}

addXStag2Bam <- function(sl.tss.index, indiv.bams.dir, stringtie.in.dir, tagXSstrandedData.path){

    premerge.bams <- list.files(
        indiv.bams.dir,
        pattern = paste0("__idx_", sl.tss.index, "\\.bam$"),
        full.names = TRUE
    )

    merged.bam <- file.path(
        stringtie.in.dir,
        paste0("idx_", sl.tss.index, "_temp.bam")
    )

    temp <- paste(
        "samtools merge -f", "-@", processors, merged.bam, paste(premerge.bams, collapse = " ")
    ) %>%
        system.cat
    
    xs.bam.file <- file.path(
        stringtie.in.dir,
        paste0("idx_", sl.tss.index, ".bam")
    )

    add.xs.tag.cmd <- paste(
        "samtools sort", "-@", processors, merged.bam, "|",
        "samtools view", "-@", processors, "-h", "|",
        " awk -v strType=1 -f ",
        file.path(tagXSstrandedData.path, "tagXSstrandedData.awk"), "|",
        "samtools view -bS -", ">",
        xs.bam.file
    )

    system.cat(add.xs.tag.cmd)

    return()
}

runStringTie <- function(sl.tss.index, stringtie.in.dir, stringtie.out.dir, all.tx.gtf){

    stringtie.out.gtf <- file.path(
        stringtie.out.dir,
        paste0("idx_", sl.tss.index, "_stringtie-out.gtf")
    )

    stringtie.cmd <- paste(
        "stringtie",
        file.path(stringtie.in.dir, paste0("idx_", sl.tss.index, ".bam")),
        "--conservative",
        "--ptf", file.path(stringtie.in.dir, paste0("promoters_idx_", sl.tss.index, ".txt")),
        "-o", stringtie.out.gtf,
        "-j", 5,
        "-p", processors,
        "-G", all.tx.gtf
    )

    system.cat(stringtie.cmd)
    return()
}


refStringTieAnnotMatching <- function(promoters.dt, all.tx.gtf, all.primary.tx.dt, stringtie.out.dir){

    ref.dt <- rtracklayer::import(all.tx.gtf) %>%
        as.data.frame %>% data.table

    ref.dt <- ref.dt[transcript_type != "retained_intron"]
    ref.dt[, exon_number := as.integer(exon_number)]
    
    readStringtieOut <- function(sl.tss.index, promoters.dt, stringtie.out.dir){
        stringtie.out.gtf <- file.path(
            stringtie.out.dir,
            paste0("idx_", sl.tss.index, "_stringtie-out.gtf")
        )

        stringtie.gr <- rtracklayer::import(stringtie.out.gtf)

        stringtie.dt <- stringtie.gr %>%
            as.data.frame %>% data.table %>% {.[, tss_index := sl.tss.index]}

        ## Fix exon_number
        stringtie.exon.dt <- stringtie.dt[type == "exon"][
            order(transcript_id, ifelse(strand == "+", 1, -1) * start)][
          , exon_number := 1:.N, by = transcript_id
        ]
        stringtie.tx.dt <- stringtie.dt[type == "transcript"]

        stringtie.dt <- rbind(
            stringtie.exon.dt,
            stringtie.tx.dt
        )[order(transcript_id, !is.na(exon_number), exon_number)]
        
        stringtie.first.exon.dt <- stringtie.exon.dt[exon_number == 1]
        
        sl.promoters.dt <- promoters.dt[tss_index == sl.tss.index] 
        sl.promoters.range.dt <- copy(sl.promoters.dt[, `:=`(
            start = start - 100,
            end = end + 100
        )])
        
        ol.dt <- findOverlaps(
            makeGRangesFromDataFrame(sl.promoters.range.dt),
            makeGRangesFromDataFrame(stringtie.first.exon.dt)
        ) %>% as.data.frame %>% data.table

        matched.tss.name2stringtie.txid.dt <- data.table(
            sl.promoters.dt[ol.dt[, queryHits], .(tss_name, q0, q25, q50, q75, q100)],
            transcript_id = stringtie.first.exon.dt[ol.dt[, subjectHits], transcript_id]
        ) 

        query.dt <- merge(
            stringtie.dt,
            matched.tss.name2stringtie.txid.dt,
            by = "transcript_id"
        ) %>%
            {.[, transcript_id := paste0(transcript_id, "_", sl.tss.index)]}

        return(query.dt)
    }

    query.dt <- mclapply(
        1:max(promoters.dt[, tss_index]),
        readStringtieOut,
        stringtie.out.dir = stringtie.out.dir,
        promoters.dt = promoters.dt,
        mc.cores = processors
    ) %>% rbindlist(use.names = TRUE, fill = TRUE)

    stringtie.score.dt <- merge(
        query.dt[type == "transcript", .(transcript_id, TPM, tss_name)],
        query.dt[
            type == "exon",
            list(coverage_len = sum(width), exon_n = max(exon_number)),
            by = list(transcript_id)
        ]
    ) %>%
        {.[, stringtie_score := as.numeric(TPM) * coverage_len / 1000]}

    ## Filter out transcripts with low score for each TSS
    stringtie.score.dt[, max_score_per_TSS := max(stringtie_score), by = tss_name]
    stringtie.score.dt <- stringtie.score.dt[stringtie_score > max_score_per_TSS * 0.5]

    query.dt <- query.dt[transcript_id %in% stringtie.score.dt[, transcript_id]]

    ## Identify the stringtie matched reference transcripts
    ref.exon.dt <- ref.dt[type == "exon"] %>%
        {.[, `:=`(
             coord_all = paste(seqnames, strand, start, end, sep = "_"),
             coord_start = paste(seqnames, strand, ifelse(strand == "+", start, end), sep = "_"),
             coord_end = paste(seqnames, strand, ifelse(strand == "+", end, start), sep = "_")
         )]}

    query.exon.dt <- query.dt[type == "exon"] %>%
        {.[, `:=`(
             coord_all = paste(seqnames, strand, start, end, sep = "_"),
             coord_start = paste(seqnames, strand, ifelse(strand == "+", start, end), sep = "_"),
             coord_end = paste(seqnames, strand, ifelse(strand == "+", end, start), sep = "_")
         )]}

    query.exon.dt[, query_all_exon_num := .N, by = transcript_id]

    setnames(
        query.exon.dt,
        old = c("exon_number", "transcript_id"),
        new = c("query_exon_number", "query_transcript_id")
    )

    ref.key.cols <- c(
        "transcript_id",
        "exon_number",
        "coord_all", "coord_start", "coord_end"
    )

    query.key.cols <- c(
        "query_transcript_id",
        "query_exon_number",
        "query_all_exon_num",
        "tss_name"
    )

    ## Identify matched reference transcripts for transcript where stringtie idengtify multiexon
    gene.body.dt <- merge(
        ref.exon.dt[, ref.key.cols, with = FALSE],
        query.exon.dt[
            query_all_exon_num != 1 &
            query_exon_number != 1 &
            query_exon_number != query_all_exon_num,
            c(query.key.cols, "coord_all"), with = FALSE],
        by = "coord_all",
        allow.cartesian = TRUE
    )

    first.exon.dt <- merge(
        ref.exon.dt[, ref.key.cols, with = FALSE],
        query.exon.dt[
            query_all_exon_num != 1 &
            query_exon_number == 1,
            c(query.key.cols, "coord_end"), with = FALSE],
        by = "coord_end"
    )

    last.exon.dt <- merge(
        ref.exon.dt[, ref.key.cols, with = FALSE],
        query.exon.dt[
            query_all_exon_num != 1 &
            query_exon_number == query_all_exon_num,
            c(query.key.cols, "coord_start"), with = FALSE],
        by = "coord_start"
    )

    all.matched.exon.dt <- rbind(
        gene.body.dt, first.exon.dt, last.exon.dt
    )


    all.matched.exon.dt[, exon_number := as.integer(exon_number)][, `:=`(
                           min_query_exon_num = min(query_exon_number),
                           ref_max_exon_num = max(exon_number),
                           ref_min_exon_num = min(exon_number),
                           ref_diff_exon_num = max(exon_number) - min(exon_number) + 1    
                       ), by = list(transcript_id, query_transcript_id)]

    ## second exon start information is useful to judge CDS canonicalness for multi_novel_first_exon
    second.exon.info.dt <- all.matched.exon.dt[
        query_exon_number == 2][
      , .(query_transcript_id, coord_all)][, `:=`(
                                   strand = str_split_fixed(coord_all, "_", n = 4)[, 2],
                                   start = str_split_fixed(coord_all, "_", n = 4)[, 3],
                                   end = str_split_fixed(coord_all, "_", n = 4)[, 4]
                               )][
      , second_exon_start := ifelse(strand == "+", start, end)
    ][,.(query_transcript_id, second_exon_start)][!duplicated(query_transcript_id)]

    all.matched.exon.dt <- merge(
        all.matched.exon.dt,
        second.exon.info.dt,
        by = "query_transcript_id"
    )

    ## Identify those where all exons from stringtie prediction matched with reference annotation
    multi.exon.matched.tx.dt <- all.matched.exon.dt[
      , list(.N, ref_diff_exon_num, query_all_exon_num, tss_name, ref_min_exon_num, second_exon_start),
        by = list(transcript_id, query_transcript_id)
    ][!duplicated(paste(transcript_id, query_transcript_id, sep = "_"))][
        N == ref_diff_exon_num & N == query_all_exon_num
    ][, exon_match := "multi_all"][
      , .(transcript_id, query_transcript_id, exon_match, tss_name, ref_min_exon_num, second_exon_start)
    ]

    ## Identify those where all exons except the first exons from stringtie prediction matched with reference annotation
    multi.novel.first.exon.matched.tx.dt <- all.matched.exon.dt[
      , list(.N, ref_diff_exon_num, query_all_exon_num, tss_name, min_query_exon_num, ref_min_exon_num, second_exon_start),
        by = list(transcript_id, query_transcript_id)
    ][!duplicated(paste(transcript_id, query_transcript_id, sep = "_"))][
        min_query_exon_num == 2 &
        (query_all_exon_num - 1) == ref_diff_exon_num &
        (query_all_exon_num - 1) == N
    ][, exon_match := "multi_novel_first_exon"][
      , .(transcript_id, query_transcript_id, exon_match, tss_name, ref_min_exon_num, second_exon_start)
    ]

    ## In case stringtie identify single exon only
    single.exon.query.dt <- copy(query.exon.dt)[
        query_all_exon_num == 1
    ] %>%
        {.[, exon_end := ifelse(strand == "+", end, start)]}

    single.exon.end.gr <- makeGRangesFromDataFrame(
        single.exon.query.dt,
        start.field = "exon_end",
        end.field = "exon_end"
    )

    ref.exon.gr <- makeGRangesFromDataFrame(ref.exon.dt)

    ol.dt <- findOverlaps(
        single.exon.end.gr,
        ref.exon.gr
    ) %>% as.data.frame %>% data.table

    single.exon.matched.tx.dt <- data.table(
        transcript_id = ref.exon.dt[ol.dt[, subjectHits], transcript_id],
        query_transcript_id = single.exon.query.dt[ol.dt[, queryHits], query_transcript_id],
        exon_match = "single",
        tss_name = single.exon.query.dt[ol.dt[, queryHits], tss_name],
        ref_min_exon_num = ref.exon.dt[ol.dt[, subjectHits], exon_number],
        second_exon_start = NA
    )

    ## Merge results for multi exons and single exon case
    exon.matched.tx.dt <- rbind(
        multi.exon.matched.tx.dt,
        multi.novel.first.exon.matched.tx.dt,
        single.exon.matched.tx.dt
    )

    ## Select transcript for each TSS
    ## Relative poistion to start codon is an important indicator
    start.codon.dt <- 
        ref.dt[
            transcript_type == "protein_coding" & type == "CDS",
            .(transcript_id, seqnames, start, end, strand)
        ] %>%
        {.[, list(
             start_codon_start = case_when(
                 strand == "+" ~ min(start),
                 strand == "-" ~ max(end)
             )
         ), by = transcript_id]} %>%
        {.[!duplicated(transcript_id)]}


    setnames(stringtie.score.dt, old = "transcript_id", new = "query_transcript_id")

    score.per.gene.dt <-  merge(
        exon.matched.tx.dt,
        stringtie.score.dt,
        by = c("query_transcript_id", "tss_name")
    )  %>% {
        merge(
            x = all.primary.tx.dt,
            y = .,
            by = "transcript_id"
        )
    } %>%
        merge(
            y = start.codon.dt,
            by = "transcript_id",
            all.x = TRUE
        )

    score.per.tss.dt <- merge(
        x = score.per.gene.dt,
        y = promoters.dt[, .(tss_name, strand, q0, q25, q50, q75, q100)],
        by = "tss_name"
    )

    score.per.tss.dt[, max_score_per_gene := max(stringtie_score), by = tss_name]
    ## compulsory filter # Start codon should be downstram of q50
    score.per.tss.dt <- score.per.tss.dt[stringtie_score > 0.8 * max_score_per_gene]

    score.per.tss.dt[
      , CDS_down_to_TSS := case_when(
            is.na(start_codon_start) ~ FALSE,
            exon_match == "multi_novel_first_exon" & strand == "+" ~
                start_codon_start >= second_exon_start,
            exon_match == "multi_novel_first_exon" & strand == "-" ~
                start_codon_start <= second_exon_start,
            strand == "+" ~ start_codon_start >= q50,
            strand == "-" ~ start_codon_start <= q50
        )]

    ## Process protein_coding and non_coding genes separately
    tx.assigned.tss.dt <- score.per.tss.dt[
        order(
            is.na(start_codon_start), #1. Prioritize protein coding tx
            !CDS_down_to_TSS, #2. start codon downstream to TSS
            factor(exon_match, levels = c("multi_all", "single", "multi_novel_first_exon")), #3. Match type
            factor(refseq_tag, levels = c("original_refseq", "refseq", "N/A")), #4. Prioritize RefSeq tx
            -(basic_tag == "basic"), #5 Prioritize basic tx
            -(utr5_len != 0),
            -(utr3_len != 0),
            -cds_len, -utr5_len, -utr3_len,
            -stringtie_score
        )][, head(.SD, 1), by = list(tss_name)]

    tx.assigned.tss.dt[
      , CDS_canonical := case_when(
            transcript_type != "protein_coding" ~ "non_coding",
            transcript_type == "protein_coding" & CDS_down_to_TSS == TRUE ~ "canonical_CDS",
            TRUE ~ "redefined_CDS"
        )]

    return(list(tx.assigned.tss.dt = tx.assigned.tss.dt, ref.dt = ref.dt, query.exon.dt = query.exon.dt))
}

createGtfForAssignedTx <- function(tx.assigned.tss.dt, ref.dt, tx.per.tss.gtf, query.exon.dt){

    ## Assign annotation to tx
    redef.ref.dt <- merge(
        tx.assigned.tss.dt[
          , .(
                tss_name, transcript_id, q0, q25, q50, q75, q100,
                ref_min_exon_num, CDS_canonical
            )],
        ref.dt,
        by = "transcript_id"
    )

    redef.ref.exon.dt <- redef.ref.dt[type == "exon"]

    ## Remove exons not identified by stringtie matching
    redef.ref.exon.dt <- redef.ref.exon.dt[exon_number >= ref_min_exon_num]
    ## Add stringtie exon for multi_novel_first_exon genes

    novel.first.exon.redef.exon.dt <- merge(
        redef.ref.exon.dt[
            tss_name %in%
            tx.assigned.tss.dt[exon_match == "multi_novel_first_exon", tss_name]
        ][exon_number == ref_min_exon_num][
          , .(tss_name, q0, q25, q50, q75, q100, ref_min_exon_num,
              source, type, score, phase,
              gene_id, gene_name, transcript_id, biotype, transcript_type,
              CDS_canonical, exon_number)
        ],
        tx.assigned.tss.dt[exon_match == "multi_novel_first_exon"][
          , .(tss_name)
        ],
        by = "tss_name"
    ) %>%
        merge(
            y = query.exon.dt[
                query_transcript_id %in% tx.assigned.tss.dt[, query_transcript_id]
            ][query_exon_number == 1][, .(tss_name, seqnames, start, end, width, strand)],
            by = "tss_name"
        )

    novel.first.exon.redef.exon.dt[, exon_number := ref_min_exon_num - 1]

    redef.ref.exon.dt <- rbind(
        redef.ref.exon.dt,
        novel.first.exon.redef.exon.dt
    )[order(tss_name, exon_number)]

    trimExonsByTSS <- function(dt){
        ## Currently use q0 for the range
        dt <- copy(dt)

        dt <- dt[
            ## Remove exons upstream of TSS
            case_when(
                strand == "+" ~ end > q0,
                strand == "-" ~ start < q100,
                TRUE ~ TRUE
            )
        ]

        dt[, first_exon_num := min(exon_number), by = tss_name]

        t.first.exon.dt <- dt[exon_number == first_exon_num & type == "exon"]
        t.first.exon.dt[, `:=`(
            ## Shift exon start position where TSS locates within the exon
            start = case_when(
                strand == "+" ~ q0,
                TRUE ~ start
            ),
            end = case_when(
                strand == "-" ~ q100,
                TRUE ~ end
            )
        )][, width := end - start + 1]
        
        dt <- rbind(
            t.first.exon.dt,
            dt[!(exon_number == first_exon_num & type == "exon")]
        )[
            order(
                tss_name,
                ifelse(strand == "+", 1, -1) * start
            )    
        ]

        dt[, first_exon_num := NULL]
        
        return(dt)
    }

    redef.ref.exon.dt <- trimExonsByTSS(redef.ref.exon.dt)

    original.gtf.cols <- colnames(redef.ref.exon.dt)

    ## swap the name for txdb transformation
    redef.ref.dt <- rbind(
        redef.ref.dt[
            type %in% c("start_codon", "stop_codon", "CDS")
        ][exon_number >= ref_min_exon_num],
        redef.ref.exon.dt 
    )[order(transcript_id, exon_number)]

    redef.ref.dt[, original_transcript_id := transcript_id][
      , `:=`(transcript_id = tss_name)
    ]

    redef.ref.exon.dt[, original_transcript_id := transcript_id][
      , `:=`(transcript_id = tss_name)
    ]

    tx.assigned.tss.dt[, original_transcript_id := transcript_id][
      , `:=`(transcript_id = tss_name)
    ]

    ## Assign annotation for tx with canonical CDS or non_coding genes
    canonical.CDS.redef.ref.dt <- redef.ref.dt[
        CDS_canonical %in% c("canonical_CDS", "non_coding")
    ] %>% trimExonsByTSS


    ## Assign annotation for tx with non-canoniccal CDS
    ## 1. Identify ORF for non-canoncal CDS transcript

    redef.ref.exon.for.non.canonical.cds.dt <- copy(redef.ref.exon.dt)[
        CDS_canonical == "redefined_CDS"
    ]

    redef.ref.exon.for.non.canonical.cds.txdb <- makeTxDbFromGRanges(
        makeGRangesFromDataFrame(
            redef.ref.exon.for.non.canonical.cds.dt,
            keep.extra.columns = TRUE
        )
    )

    library("BSgenome.Hsapiens.UCSC.hg38")

    ref.exon.grl <- exonsBy(redef.ref.exon.for.non.canonical.cds.txdb, by = "tx")

    names(ref.exon.grl) <- AnnotationDbi::select(
                                              redef.ref.exon.for.non.canonical.cds.txdb,
                                              keys = names(ref.exon.grl),
                                              columns = "TXNAME", keytype = "TXID"
                                          )[["TXNAME"]]


    redef.ref.exon.for.non.canonical.cds.seq <- extractTranscriptSeqs(
        BSgenome.Hsapiens.UCSC.hg38,
        ref.exon.grl
    )

    library("ORFik")

    ref.orf.dt <- findMapORFs(
        ref.exon.grl,
        redef.ref.exon.for.non.canonical.cds.seq,
        startCodon = "ATG",
        longestORF = TRUE,
        minimumLength = 90,
        groupByTx = TRUE
    ) %>% as.data.frame %>% data.table

    setnames(ref.orf.dt, old = c("group_name", "names"), new = c("tss_name", "orf_name"))

    extract.ref.orf.stop.dt <- ref.orf.dt[, list(
        tss_name,
        start_codon_start = case_when(
            strand == "+" ~ min(start),
            strand == "-" ~ max(end)
        ),    
        stop_codon_end = case_when(
            strand == "+" ~ max(end),
            strand == "-" ~ min(start)
        ),
        CDS_width = sum(width)
    ), by = "orf_name"
    ][!duplicated(orf_name)]

    ref.stop.codon.dt <- 
        redef.ref.dt[
            transcript_type == "protein_coding" & type == "stop_codon",
            .(transcript_id, seqnames, start, end, strand)
        ] %>%
        {.[, list(
             stop_codon_end = case_when(
                 strand == "+" ~ max(end),
                 strand == "-" ~ min(start)
             )
         ), by = transcript_id]} %>%
        {.[!duplicated(transcript_id)]}

    assigned.orf.dt <- merge(
        tx.assigned.tss.dt[
            CDS_canonical == "redefined_CDS",
            .(tss_name, transcript_id, strand, q0, q25, q50, q75, q100)],
        ref.stop.codon.dt,
        by = "transcript_id" 
    ) %>%
        merge(
            y = extract.ref.orf.stop.dt,
            by = c("tss_name", "stop_codon_end")
        )

    assigned.orf.coord.dt <- merge(
        ref.orf.dt,
        assigned.orf.dt[, .(orf_name, transcript_id, start_codon_start)],
        by = "orf_name"
    )

    assigned.orf.coord.dt[
      , `:=`(
            orf_name = NULL,
            group = NULL,
            type = "CDS",
            width = end - start + 1
        )]

    ## add exon_number information
    orf.ol.dt <- findOverlaps(
        makeGRangesFromDataFrame(redef.ref.exon.for.non.canonical.cds.dt),
        makeGRangesFromDataFrame(assigned.orf.coord.dt)
    ) %>% as.data.frame %>% data.table

    to.match.dt <- redef.ref.exon.for.non.canonical.cds.dt[
        orf.ol.dt[, queryHits],
        .(tss_name, q0, q25, q50, q75, q100, ref_min_exon_num, CDS_canonical, source, score, phase,
          gene_id, gene_name, biotype, transcript_type, exon_number, original_transcript_id)
    ]

    setnames(to.match.dt, old = "tss_name", new = "match_tss_name")

    assigned.orf.coord.dt <- cbind(
        assigned.orf.coord.dt[orf.ol.dt[, subjectHits]],
        to.match.dt
    )[tss_name == match_tss_name][, match_tss_name := NULL]

    ## Merge exon and CDS information
    noncanonical.CDS.redef.ref.dt <- merge(
        redef.ref.dt[CDS_canonical == "redefined_CDS"],
        assigned.orf.dt[, .(transcript_id, start_codon_start)],
        by = "transcript_id"
    )

    noncanonical.CDS.redef.noCDS.ref.dt <- trimExonsByTSS(noncanonical.CDS.redef.ref.dt)[
        type %in% c("exon", "stop_codon") &
        tss_name %in% unique(assigned.orf.coord.dt[, tss_name])
    ]

    noncanonical.CDS.redef.ref.dt <- rbind(
        noncanonical.CDS.redef.noCDS.ref.dt,
        assigned.orf.coord.dt,
        use.names = TRUE
    )[order(transcript_id, exon_number)]

    addStartCodon <- function(sl.tss.name, noncanonical.CDS.redef.ref.dt){
        ## print(sl.tss.name)
        sl.dt <- noncanonical.CDS.redef.ref.dt[tss_name == sl.tss.name]

        sl.cds.dt <- sl.dt[order(ifelse(strand == "+", 1, -1) * start)][type == "CDS"]

        sl.cds.dt[, `:=`(
            CDS_id = 1:.N
        )]

        if(sl.cds.dt[CDS_id == 1, width] > 2){
            start.codon.dt <- copy(sl.cds.dt[CDS_id == 1])
            start.codon.dt[, `:=`(
                start = ifelse(strand == "+", start, end - 2),
                end = ifelse(strand == "+", start + 2, end),
                type = "start_codon",
                width = 3
            )]
        } else {
            start.codon.1.dt <- copy(sl.cds.dt[CDS_id == 1])
            start.codon.2.dt <- copy(sl.cds.dt[CDS_id == 2])
            start.codon.2.dt[, `:=`(
                start = ifelse(strand == "+", start, end - 2 + start.codon.1.dt[, width]),
                end = ifelse(strand == "+", start + 2 - start.codon.1.dt[, width], end),
                width = 3 - start.codon.1.dt[, width]
            )]
            start.codon.dt <- rbind(start.codon.1.dt, start.codon.2.dt)
            start.codon.dt[, type := "start_codon"]
        }

        start.codon.dt[, `:=`(CDS_id = NULL)]

        sl.redefined.dt <- rbindlist(list(
            sl.dt[type != "CDS"],
            sl.cds.dt[, `:=`(CDS_id = NULL)],
            start.codon.dt
        ))[
            order(
                exon_number,
                factor(type, levels = c("exon", "start_codon", "stop_codon", "CDS"))
            )]

        return(sl.redefined.dt)
    }

    noncanonical.start.CDS.redef.ref.dt <- mclapply(
        unique(noncanonical.CDS.redef.ref.dt[type == "CDS", tss_name]),
        addStartCodon,
        noncanonical.CDS.redef.ref.dt = noncanonical.CDS.redef.ref.dt,
        mc.cores = processors
    ) %>% rbindlist

    tss.redef.ref.dt <- rbind(
        canonical.CDS.redef.ref.dt[
            CDS_canonical == "canonical_CDS",
            c(original.gtf.cols, "tss_name", "CDS_canonical"), with = FALSE],
        canonical.CDS.redef.ref.dt[
            CDS_canonical == "non_coding" & type == "exon",
            c(original.gtf.cols, "tss_name", "CDS_canonical"), with = FALSE],
        noncanonical.start.CDS.redef.ref.dt[
          , c(original.gtf.cols, "tss_name", "CDS_canonical"), with = FALSE]
    )

    ## Add phase information
    tss.redef.ref.cds.dt <- tss.redef.ref.dt[type == "CDS"]

    tss.redef.ref.cds.dt <- tss.redef.ref.cds.dt[order(transcript_id, exon_number)]
    tss.redef.ref.cds.dt[
      , cum_sum_m1 := cumsum(shift(width, fill = 0, type = "lag")), by = transcript_id
    ][, phase := mod(3 - mod(cum_sum_m1, 3), 3)][, cum_sum_m1 := NULL]

    tss.redef.ref.with.phase.dt <- rbind(
        tss.redef.ref.dt[type != "CDS"],
        tss.redef.ref.cds.dt
    )[order(seqnames, transcript_id, exon_number)]


    tss.redef.ref.gr <- makeGRangesFromDataFrame(
        tss.redef.ref.with.phase.dt,
        keep.extra.columns = TRUE
    ) %>% sort

    rtracklayer::export(
                     tss.redef.ref.gr,
                     tx.per.tss.gtf
                 )

    return(tss.redef.ref.dt)
}
