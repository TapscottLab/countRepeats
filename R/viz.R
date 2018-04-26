.getSigRepName <- function(ddsRepName, padj.thres=0.1, lfc.thres=0) {
    res <- results(ddsRepName)
    res <- cbind(rowData(ddsRepName), res)
    keep <- which(res$padj < padj.thres & abs(res$log2FoldChange) > lfc.thres)
    sig <- res[keep, ]
}
#'
#' visualization of two groups of log values (log2(normalized count + 0.5))
#'
scatterLog2PerGroup <- function(ddsRepName, factor,
                                padj.thres=0.1, lfc.thres=0,
                                output.path=".") {
    #' sanity check
    if (!is.factor(factor))
        stop("factor argument must be a factor instance")
    lev <- levels(factor)
    if (length(lev) > 2)
        stop("Scatter plot only allows two groups of samples. Check factor.")
    if (length(lev) == 1)
        stop("Scatter plot need to have two groups of samples. Check factor")

    res <- results(ddsRepName)
    sig <- res$padj < padj.thres &  abs(res$log2FoldChange) > lfc.thres
    sig[is.na(sig)] <- FALSE
    cnt <- counts(ddsRepName, normalized=TRUE)
    log2.cnt <- log2(cnt+.5)
    df <- data.frame(y=rowMeans(log2.cnt[, factor==lev[2]]),
                     x=rowMeans(log2.cnt[, factor==lev[1]]),
                     sig=factor(sig))
    mycol <- c("steelblue", "red")
    title <- paste0("rmskStats_Log2ScatterPlot_", lev[2], "vs", lev[1])
    ggfile <- file.path(output.path, paste0(title, ".pdf"))
    gg <- ggplot(df, aes(x=x, y=y), alpha=0.5) +
        geom_point(col=mycol[df$sig], alpha=0.5) +
          theme_classic() + geom_smooth(method="lm", col="grey") +
            labs(x=lev[1], y=lev[2], title="log2(normalized count + 0.5)")
    ggsave(file=ggfile, gg)

}

#'
#' Heatmap of repName
#'
heatmapRepName <- function(ddsRepName, factor,
                           sig.only=TRUE,
                           padj.thres=0.1, lfc.thres=0,
                           title=NULL, row.annotation=TRUE,
                           output.path=".", ...) {
    require(pheatmap)
    file <- file.path(output.path, paste0("rmskStats_Heatmap_repName_", title, ".pdf"))
    rlg <- rlog(ddsRepName)
    if (sig.only) {
        targets <- .getSigRepName(ddsRepName, padj.thres=padj.thres,
                                  lfc.thres=lfc.thres)
        rlg <- rlg[rownames(targets)]
    }
    mat <- assay(rlg)

    #' scale by center
    annotation_col <- data.frame(group=factor)
    rownames(annotation_col) <- colnames(mat)
    annotation_row <- data.frame(Family=rowData(rlg)$repFamily,
                                 Class=rowData(rlg)$repClass)
    rownames(annotation_row) <- rownames(mat)

    pheatmap(mat, #annotation_row=annotation_row,
             #annotation_col=annotation_col,
             scale="row",
             silent=TRUE,
             filename=file, ...)

}
    
#'
#' visualize the number of significent repNames
#'

vizSigFreqOnFamily <- function(ddsRepName, padj.thres=0.1, lfc.thres=0.5, main=NULL, output.path=".") {
    require(DESeq2)
    res <- results(ddsRepName)
    res <- cbind(rowData(ddsRepName), res)
    keep <- which(res$padj < padj.thres & abs(res$log2FoldChange) > lfc.thres)
    sig <- res[keep, ]
    sig$expression <- ifelse(sig$log2FoldChange > 0, "up", "down")
    sig$repFamily <- factor(sig$repFamily)
    sig$repClass <- factor(sig$repClass)
    sig$nick_name <- main
    sig <- as.data.frame(sig)
    title <- sprintf("%s  %3.0f and %3.0f repNames up- and down-regulated",
                     main, sum(sig$expression == "up"),
                     sum(sig$expression== "down"))
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_SigRepName_on_repFamily.pdf"))
    gg <- ggplot(sig, aes(x=factor(repFamily), fill=expression)) +
      geom_bar(stat="count", width=0.7) +
        labs(title=title, x="repFamily", y="frequency of sig. repNames") +
            theme(axis.text.x=element_text(angle=90, hjust=1))
    ggsave(file=ggfile, gg)
    return(gg)
}

vizSigFreqOnClass <- function(ddsRepName, padj.thres=0.1, fc.thres=0.5, main=NULL, output.path=".") {
    require(DESeq2)
    res <- results(ddsRepName)
    res <- cbind(rowData(ddsRepName), res)
    keep <- which(res$padj < padj.thres & abs(res$log2FoldChange) > fc.thres)
    sig <- res[keep, ]
    sig$expression <- ifelse(sig$log2FoldChange > 0, "up", "down")
    sig$repFamily <- factor(sig$repFamily)
    sig$repClass <- factor(sig$repClass)
    sig$nick_name <- main
    sig <- as.data.frame(sig)
    title <- sprintf("%s  %3.0f and %3.0f repNames up- and down-regulated",
                     main, sum(sig$expression == "up"),
                     sum(sig$expression== "down"))
    gg <- ggplot(sig, aes(x=factor(repClass), fill=expression)) +
      geom_bar(stat="count", width=0.7) +
        labs(title=title, x="repClass", y="frequency of sig. repNames") +
            theme(axis.text.x=element_text(angle=90, hjust=1))
    
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_SigRepName_on_repClass.pdf"))
    ggsave(file=ggfile, gg)
    return(gg)
}


vizEnriched<- function(rmsk_stats, main=NULL, rmsk.pval=0.1, output.path=".") {
    ## repstats <- rmsk.res[[x]]$repstats
    ## family
    enriched_fam <- as.data.frame(rmsk_stats$fam_enrichment)
    enriched_fam <- subset(enriched_fam, num.sigRepName > 0)
    enriched_fam$repFamily <- rownames(enriched_fam)
    enriched_fam$enrichment <- ifelse(enriched_fam$prob < rmsk.pval &
                                          enriched_fam$num.sigRepName > enriched_fam$mu,
                                      "Yes", "No")
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_enriched_repFam.pdf"))
    gg1 <- ggplot(enriched_fam, aes(x=repFamily, y=num.sigRepName, fill=enrichment)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repFamily, y=mu), color="red") +
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=main, y="frequency of sig. repName") +
                    scale_fill_manual(values=c("#999999", "#56B4E9"))
    ggsave(file=ggfile, gg1)

    ## class
    enriched_class <- as.data.frame(rmsk_stats$class_enrichment)
    enriched_class <- subset(enriched_class, num.sigRepName > 0)
    enriched_class$repClass <- rownames(enriched_class)
    enriched_class$enrichment <- ifelse(enriched_class$prob < rmsk.pval &
                                          enriched_class$num.sigRepName > enriched_class$mu,
                                      "Yes", "No")
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_enriched_repClass.pdf"))
    gg2 <- ggplot(enriched_class, aes(x=repClass, y=num.sigRepName, fill=enrichment)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repClass, y=mu), color="red") +
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=main, y="frequency of sig. repName") +
                    scale_fill_manual(values=c("#999999", "#56B4E9"))
    ggsave(file=ggfile, gg2)
    return(list(fam=gg1, class=gg2))
    
}
##vizDepleted

vizDepleted <- function(rmsk_stats, main=NULL, rmsk.pval=0.1, output.path=".") {
    ## repstats <- rmsk.res[[x]]$repstats
    ## family
    depleted_fam <- as.data.frame(rmsk_stats$fam_depletion)
    depleted_fam <- subset(depleted_fam, num.sigRepName > 0)
    depleted_fam$repFamily <- rownames(depleted_fam)
    depleted_fam$depletion <- ifelse(depleted_fam$prob < rmsk.pval &
                                          depleted_fam$num.sigRepName > depleted_fam$mu,
                                      "Yes", "No")
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_depleted_repFam.pdf"))
    gg1 <- ggplot(depleted_fam, aes(x=repFamily, y=num.sigRepName, fill=depletion)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repFamily, y=mu), color="red") +
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=main, y="frequency of sig. repName") +
                    scale_fill_manual(values=c("#999999", "#56B4E9"))
    ggsave(file=ggfile, gg1)

    ## class
    depleted_class <- as.data.frame(rmsk_stats$class_depletion)
    depleted_class <- subset(depleted_class, num.sigRepName > 0)
    depleted_class$repClass <- rownames(depleted_class)
    depleted_class$depletion <- ifelse(depleted_class$prob < rmsk.pval &
                                        depleted_class$num.sigRepName > depleted_class$mu,
                                      "Yes", "No")
    ggfile <- file.path(output.path, paste0("rmsk_", main, "_depleted_repClass.pdf"))
    gg2 <- ggplot(depleted_class, aes(x=repClass, y=num.sigRepName, fill=depletion)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repClass, y=mu), color="red") +
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=main, y="frequency of sig. repName") +
                    scale_fill_manual(values=c("#999999", "#56B4E9"))
    ggsave(file=ggfile, gg2)
    return(list(fam=gg1, class=gg2))
    
}
##vizDepleted
