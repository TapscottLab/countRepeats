summarizeRepNameCount <- function(rep_bams, rmsk, paired_end=TRUE, cores=2L) {
    require(rmskStats)
    require(GenomicAlignments)
    require(parallel)
    rmsk_grl <- rmskStats::splitByRepName(rmsk)

    cnts <- mclapply(rep_bams, function(x) {
        message(basename(x))
        flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
        what <- c("qname", "flag", "qwidth", "isize", "mapq")
        param <- ScanBamParam(flag = flag, what = what, tag = c("NM", "XA", "XS"))
        if (paired_end)
            ga <- readGAlignmentPairs(x, param=param)
        else
            ga <- readGAlignments(x, param=param)
      
        tb <- table(seqnames(ga))
    }, mc.cores=cores, mc.preschedule=FALSE)
    cnts <- do.call(cbind, cnts)
    colnames(cnts) <- sub(".bam", "", basename(rep_bams))

    #' sanity check
    features <- rownames(cnts)[rownames(cnts) %in% names(rmsk_grl)]
    cnts <- cnts[features, ]
    rmsk_grl <- rmsk_grl[features, ]
    #' make SummarizedExperiment instance
    se <- SummarizedExperiment(assays=SimpleList(counts=cnts),
                               rowRanges=rmsk_grl)
    se
}  
