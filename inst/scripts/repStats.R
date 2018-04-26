.getClassFamily <- function(gr) {
        c(repClass=paste(as.character(unique(gr$repClass)), collapse=";"),
          repFamily=paste(as.character(unique(gr$repFamily)), collapse=";"))
}

getRepNameFeatures <- function(rmsk) {
    features <- split(rmsk, rmsk$repName)

    #' fix (CATTC)n and (GAATC)n - split into two sattellite, simple_repeat
    tmp <- split(features[["(CATTC)n"]], factor(features[["(CATTC)n"]]$repClass))
    features[["(CATTC)n-Satellite"]] <- tmp[["Satellite"]]
    features[["(CATTC)n-Simple_repeat"]] <- tmp[["Simple_repeat"]]

    tmp <- split(features[["(GAATG)n"]], features[["(GAATG)n"]]$repClass)
    features[["(GAATG)n-Satellite"]] <- tmp[["Satellite"]]
    features[["(GAATG)n-Simple_repeat"]] <- tmp[["Simple_repeat"]]

    i <- which(names(features) %in% c("(CATTC)n", "(GAATG)n"))
    features <- features[-i]
    return(features)
}

countRepeatName <- function(rmsk, all.reads, cores=1) {
    features <- split(rmsk, rmsk$repName)
    ## Some repName are associated with two repFamily, so split them up
    ## (CATTC)n and (GAATG)n
    repNameSE <- summarizeOverlaps.adjNH(features=features,
                                         all.reads=all.reads,
                                         type="within",
                                         inter.feature=TRUE,
                                         ignore.strand=TRUE,
                                         cores=cores)
}

repDESeq2 <- function(se, formula, ...) {
    require(DESeq2)
    dds <- DESeqDataSet(se, design = formula)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds, ...)
    dds 
}

repNameSig <- function(dds.repName, padj.thres=0.05, fc.thres=0.95) {
    ## dds.repName: DESeq2 dataset
    res <- results(dds.repName)
    sig <- which(res$padj < padj.thres & abs(res$log2FoldChange) > fc.thres)
    dds.repName[sig]
}

do.hyper <- function(universe, selected, x, k) {
    m <- selected
    n <- universe - selected
    prob <- dhyper(x=x, m=m, n=n, k=k)
    mu <- as.numeric(k * (m/(m+n)))
    DF <- DataFrame(prob=prob,
                    mu=mu,
                    num.sigRepName=as.integer(x),
                    num.totalRepName=as.integer(k))
    return(DF)
}

repStats <- function(se, sigRepName) {
    sigResults <- results(sigRepName)
    sigEnriched <- sigRepName[sigResults$log2FoldChange > 0, ]
    sigDepleted <- sigRepName[sigResults$log2FoldChange < 0, ]

    universe <- length(se)
    sel_enriched <- length(sigEnriched)
    sel_depleted <- length(sigDepleted)
    
    ## family enrichement
    k <- tb.all <- table(factor(mcols(se)$repFamily))
    tb.sig <- table(factor(mcols(sigEnriched)$repFamily))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    fam_enrichment <- do.hyper(universe, selected=sel_enriched,
                                x=x, k=k)

    ## family depletion
    tb.sig <- table(factor(mcols(sigDepleted)$repFamily))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    fam_depletion <- do.hyper(universe, selected=sel_depleted,
                              x=x, k=k)

    ## Class enrichemnt
    k <- table(factor(mcols(se)$repClass))
    tb.sig <- table(factor(mcols(sigEnriched)$repClass))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    class_enrichment <- do.hyper(universe, selected=sel_enriched,
                                 x=x, k=k)
    ## Class depletion
    tb.sig <- table(factor(mcols(sigDepleted)$repClass))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    class_depletion <- do.hyper(universe, selected=sel_depleted,
                               x=x, k=k)

    return(list(fam_enrichment=fam_enrichment,
                fam_depletion=fam_depletion,
                class_enrichment=class_enrichment,
                class_depletion=class_depletion))
}

summarizeRepStats <- function(repstats, pval=0.05, verbose=TRUE) {
    res <- lapply(repstats, function(x) x[x$prob < pval & x$num.sig > x$mu, ,
                                          drop=FALSE]) 
    if (verbose) {
        message("Enriched Family:")
        print(res[[1]])
        message("\n\nDepleted Family:")
        print(res[[2]])
        message("\n\nEnriched Class:")
        print(res[[3]])
        message("\n\nDepleted Class:")
        print(res[[4]])
    }
    return(res)
}

repStatsReport <- function(ddsRepName, sigRepName, repstats, file,
                           pval=0.05) {
    require(xlsx)
    append <- FALSE
    res <- results(ddsRepName)
    ## sheet 1: all repeat names
    out <- cbind(mcols(ddsRepName)[, c("repFamily", "repClass")],
                 res,
                 as(assay(ddsRepName), "DataFrame"))
    
    if (nrow(out) > 65536) { ## row limit for worksheet via Java
        message("Put the whole repNames statistics on ddsRepName_statisitcs.csv")
        f <- file.path(dirname(file), "ddsRepName_statistics.csv")
        write.csv(out, file=f)
    }
    else { 
        write.xlsx(as.data.frame(out), file=file, sheetName="All RepeatName",
                   append=FALSE)
        append <- TRUE
    }

    ## sheet 2: sig repeat names
    sig_dds <- ddsRepName[rownames(sigRepName)]
    sig_res <- res[rownames(sigRepName), ]
    out <- cbind(mcols(sig_dds)[, c("repFamily", "repClass")],
                 sig_res,
                 as(assay(sigRepName), "DataFrame"))
    write.xlsx(as.data.frame(out), file=file, sheetName="Sig RepeatName",
               append=append)
    
    ## sheet 3: result of Family/class enrichment/depletion analysis
    res <- summarizeRepStats(repstats, verbose=FALSE, pval=pval)
    start_row <- 1
    sheetname <- "RFEA Results"
    wb <- loadWorkbook(file)
    sheets <- getSheets(wb)
    
    if (any(names(sheets) %in% sheetname)) removeSheet(wb, sheetname)
    
    sheet <- createSheet(wb, sheetName=sheetname)
    for (i in names(res)) {
        addDataFrame(as.data.frame(i), sheet=sheet,
                     col.names=FALSE, row.names=FALSE,
                     startRow=start_row)
        start_row <- start_row + 1
        if (nrow(res[[i]]) > 0) {
        addDataFrame(as.data.frame(res[[i]]), sheet=sheet,
                     startRow=start_row)
        start_row <- start_row + nrow(res[[i]])+3
        }
    }

    saveWorkbook(wb, file=file)

    ## sheet 4 to 7
    lapply(names(repstats), function(x, file) {
        if (nrow(repstats[[x]]) > 0)
        write.xlsx(as.data.frame(repstats[[x]]),
                                 file=file, sheetName=x,
                                 append=TRUE)
    }, file)
    return(invisible())
}

do.repStats <- function(se, file) {
    formula <- as.formula("~ Treatment")
    ddsRepName <- repDESeq2(se, formula)
    sigRepName <- repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
    print(sigRepName)
    repstats <- repStats(ddsRepName, sigRepName)
    res <- summarizeRepStats(repstats)
    repStatsReport(ddsRepName, sigRepName, repstats,
                   file=file, pval=0.05)
    res
}

do.repStats <- function(se, output.file, formula) {
    #formula <- as.formula("~ Treatment")
    ddsRepName <- repDESeq2(se, formula)
    sigRepName <- repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
    print(sigRepName)
    repstats <- repStats(ddsRepName, sigRepName)
    res <- summarizeRepStats(repstats)
    repStatsReport(ddsRepName, sigRepName, repstats,
                   file=output.file, pval=0.05)
    res
}
