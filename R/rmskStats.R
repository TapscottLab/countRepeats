.do_hyperGeom <- function(universe, selected, x, k) {
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

rmskStats <- function(dds, alpha=0.1, lfc.threshold=NULL, verbose=TRUE) {
    res <- results(dds, alpha=alpha)
    sigResults <- subset(res, res$padj < alpha )
    
    if (!is.null(lfc.threshold))
        sigResults <- subset(sigResults, abs(sigResults$log2FoldChange) > lfc.threshold)
    
    up   <- rownames(sigResults)[sigResults$log2FoldChange > 0]
    down <- rownames(sigResults)[sigResults$log2FoldChange < 0]

    if (verbose) {
        msg <- sprintf("Total %g repNames are siganificant with lfc > %0.2f and padj < %0.2f",
                       nrow(sigResults), lfc.threshold, alpha)
        message(msg)
    }
    ## split sig repNames by enriched and depleted
    sigEnriched <- dds[up, ]
    sigDepleted <- dds[down, ]

    universe <- length(dds)
    sel_enriched <- length(sigEnriched)
    sel_depleted <- length(sigDepleted)
    
    ## family enrichement
    k <- tb.all <- table(factor(mcols(dds)$repFamily))
    tb.sig <- table(factor(mcols(sigEnriched)$repFamily))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    fam_enrichment <- .do_hyperGeom(universe, selected=sel_enriched,
                                x=x, k=k)

    ## family depletion
    tb.sig <- table(factor(mcols(sigDepleted)$repFamily))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    fam_depletion <- .do_hyperGeom(universe, selected=sel_depleted,
                              x=x, k=k)

    ## Class enrichemnt
    k <- table(factor(mcols(dds)$repClass))
    tb.sig <- table(factor(mcols(sigEnriched)$repClass))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    class_enrichment <- .do_hyperGeom(universe, selected=sel_enriched,
                                 x=x, k=k)
    ## Class depletion
    tb.sig <- table(factor(mcols(sigDepleted)$repClass))
    x <- rep(0, length(k))
    names(x) <- names(k)
    x[names(tb.sig)] <- tb.sig
    class_depletion <- .do_hyperGeom(universe, selected=sel_depleted,
                               x=x, k=k)
    #' should result an instance called rmskStatsResult with log: alpha and lfc.threshold
    return(list(fam_enrichment=fam_enrichment,
                fam_depletion=fam_depletion,
                class_enrichment=class_enrichment,
                class_depletion=class_depletion,
                alpha=alpha,
                lfc.threshold=lfc.threshold))
}

.summary <- function(rmsk_stats, rmsk.pval=0.1, verbose=TRUE) {
    #' dfl is a list of DataFrame: should change to "rmskStatsRestult"
    res <- lapply(rmsk_stats[1:4], function(x)
        subset(x, x$prob < rmsk.pval & x$num.sigRepName > x$mu))
    cat("Enrichment/Depletion P-value threshold:", rmsk.pval,
        "\n")
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

rmskStatsReport <- function(ddsRepName, rmsk_stats, file,
                            rmsk.pval=0.1, verbose=FALSE) {
    require(xlsx)
    res <- results(ddsRepName, alpha=rmsk_stats$alpha)
    append <- FALSE
    sig_res <- subset(res, res$padj < rmsk_stats$alpha)
    
    if (!is.null(rmsk_stats$lfc.threshold))
        sig_res <- subset(sig_res, abs(sig_res$log2FoldChange) > rmsk_stats$lfc.threshold)
    
    sig_dds <- ddsRepName[rownames(sig_res)]
    
    ## sheet 1: all repeat names
    message("Sheet 1: all repeat names")
    norm.count <- counts(ddsRepName, normalized=TRUE)
    colnames(norm.count) <- paste0("norm_", colnames(norm.count))
    out <- cbind(mcols(ddsRepName)[, c("repFamily", "repClass")],
                 res,
                 as(counts(ddsRepName), "DataFrame"),
                 as(norm.count, "DataFrame"))
    out <- out[order(out$padj, decreasing=FALSE), ]
    if (nrow(out) > 65536) { ## row limit for worksheet via Java
        message("Put the whole repNames statistics on ddsRepName_statisitcs.csv")
        f <- file.path(dirname(file), "ddsRepName_statistics.csv")
        write.csv(out, file=f)
    }
    else { 
        write.xlsx(as.data.frame(out), file=file, sheetName="RepName DESeq Results ",
                   append=FALSE)
        append <- TRUE
    }

    ## sheet 2: sig repeat names
    message("Sheet 2: significant repeat names")
    out <- cbind(mcols(sig_dds)[, c("repFamily", "repClass")],
                 sig_res,
                 as(counts(sig_dds), "DataFrame"))
    out <- out[order(out$padj, decreasing=FALSE), ]
    write.xlsx(as.data.frame(out), file=file, sheetName="Sig RepeatName",
               append=append)
    
    ## sheet 3: result of Family/class enrichment/depletion analysis
    res <- .summary(rmsk_stats, verbose=FALSE, rmsk.pval=rmsk.pval)
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
    lapply(names(rmsk_stats[1:4]), function(x, file) {
        message(x)
        if (nrow(rmsk_stats[[x]]) > 0)
        write.xlsx(as.data.frame(rmsk_stats[[x]]),
                                 file=file, sheetName=x,
                                 append=TRUE)
    }, file)
    return(invisible())
}

rmskStatsWrapper <- function(se, formula, alpha=0.1, lfc.threshold=NULL,
                             file=NULL, rmsk.pval=0.1, size.factor=NULL) {
    #' (1) DESeq
    #se <- se[rowSums(assays(se)[["counts"]]) > ncol(se), ]
    dds <- DESeqDataSet(se, design = formula)
    
    if (!is.null(size.factor)) ## need sanity check in the feature
        sizeFactors(dds) <- size.factor
    
    dds <- DESeq(dds)
    rmsk_stats <- rmskStats(dds, alpha=alpha, lfc.threshold=lfc.threshold)
    #' summary - show on the workspace
    smry <- .summary(rmsk_stats, rmsk.pval=rmsk.pval, verbose=TRUE)
    #' report:
    message("Generating report:", file)
    if (!is.null(file)) {
        #' sanity check
        if (!file.exists(dirname(file)))
            message(dirname(file), "does not exist. Report file cannot be printed")
        else    
            rmskStatsReport(ddsRepName=dds, rmsk_stats=rmsk_stats,
                            file=file, rmsk.pval=rmsk.pval,
                            verbose=FALSE)
    }
    #' visualization *.pdf: scatter plots, heatmap plots

    return(list(ddsRepName=dds, rmsk_stats=rmsk_stats))
}


