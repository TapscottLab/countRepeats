#'
#' adapted from GenomicAlignment::summarizeOverlaps. Same approach but
#' the counts are adjusted by NH: number of reported alignment.
#' Only support mode: IntersectionStrict and inter.feature:TRUE/FALSE
#' 

.countSubjectHits <- function(ov, features, reads, round.up=TRUE) {
    #' if counting for larger element such as repName
    #' return integer
    NH <- mcols(reads)$NH
    fracReadCount <- 1/NH[queryHits(ov)]
    cnt <- rep(0, length(features))
    names(cnt) <- as.character(1:length(features))
    tmp <- tapply(fracReadCount, subjectHits(ov), sum)
    cnt[as.integer(names(tmp))] <- tmp
    if (round.up)  {
        cnt <- as.integer(round(cnt))
    }
    return(cnt)
}

countHits <- function(features, reads, type=c("any", "start", "end", "within"),
                  ignore.strand=TRUE, inter.feature=FALSE, round.up=TRUE)
{   #' at least overlap with 20L
    ov <- findOverlaps(reads, features,
                       type=match.arg(type),
                       minoverlap=20L,
                       ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    .countSubjectHits(ov, features=features, reads=reads,
                      round.up=round.up)
}

.bamToGAlignments <- function(bam_file, singleEnd=TRUE, strandMode=1,
                              ignore.strand=TRUE) {
    #' take the bam_file and return GAlignmentPairs or GAlignments instance
    flag <- scanBamFlag(isUnmappedQuery = FALSE,
                        isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth")
    param <- ScanBamParam(flag = flag, what = what, tag = c("XS", "NH"))

    if (singleEnd) {
        ga <- readGAlignments(bam_file, param=param)
    }
    
    if (!singleEnd) {
        #' updating flag and param
        flag2 <- scanBamFlag(isProperPair=TRUE)
        param <- ScanBamParam(flag = bamFlagAND(flag, flag2),
                              what = what, tag = c("XS", "NH"))
        #' give mcols to ga with first reads' NH
        ga <- readGAlignmentPairs(bam_file, param=param,
                                  strandMode=strandMode)
        mcols(ga)$NH <- mcols(first(ga))$NH
    }
    return(ga)
}

.processFeatures <- function(features, ignore.strand=TRUE) {
    if (ignore.strand) {
        if (class(features) == "GRangesList") {
            r <- unlist(features)
            strand(r) <- "*"
            features@unlistData <- r
        } else {
            strand(features) <- "*"
        }
    }
    return(features)
}

countRepeats <- function(features, bam_files, cores=1,
                        singleEnd=TRUE, strandMode=1,
                        type=c("any", "start", "end", "within"),
                        ignore.strand=TRUE,
                        inter.feature=FALSE,
                        round.up=TRUE) {
    ## This function is adapted from GenomicAlignments::summarizeOverlaps,
    ## with a few changes: the input must be bam files and the counts are
    ## adjusted by NH (number of reported multiple alignments).
    ##
    ## features: GRangeList or GRanges
    ## reads: a vector of character or BamFileList or BamFile
    ## type: any -> union, within -> IntersectStrict
    ## strandMode = 0, 1, 2
    type <- match.arg(type)

    if (strandMode == 0)
        ignore.strand <- TRUE

    features <- .processFeatures(features, ignore.strand=ignore.strand)
    
    counts <- BiocParallel::bplapply(bam_files, function(bam_file) {
        message(basename(bam_file))
        ga <- .bamToGAlignments(bam_file, strandMode=strandMode,
                                ignore.strand=ignore.strand)
        countHits(features=features, reads=ga,
                  type=type,
                  ignore.strand=ignore.strand,
                  inter.feature=inter.feature,
                  round.up=round.up)
    },  BPPARAM=MulticoreParam(worker=cores))

    counts <- as.matrix(do.call(cbind, counts))
    se <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowRanges=features)
    colData <- DataFrame(sample_name = sub(".bam", "", basename(bam_files)),
                         file_bam = bam_files)
    colnames(se) <- rownames(colData) <- colData$sample_name
    colData(se) <- colData
    se
}


