.finerSplitByRepeatFamily <- function(features) {
    family <- sapply(features, function(gr) unique(gr$repFamily))
    i <- which(elementNROWS(family) > 1)
    if (!identical(i, integer(0))) {
        #' split by repFamily
        tmp <- unlist(features[i], use.names=FALSE)
        new_name <- paste(tmp$repName, tmp$repFamily, sep="-")
        splited_features <- S4Vectors::split(tmp, factor(new_name))
        features <- features[-i]
        features <- c(features, splited_features)
    }
    return(features)
}

.getElementMetaData <- function(grl) {
    meta_data <- sapply(grl, function(gr) c(as.character(gr$repFamily)[1],
                                            as.character(gr$repClass)[1]))
    meta_data <- as(t(meta_data), "DataFrame")
    names(meta_data) <- c("repFamily", "repClass")
    meta_data

}

splitByRepName <- function(rmsk, finerSplitByRepFamily=TRUE) {
    features <- S4Vectors::split(rmsk, rmsk$repName)
    
    #' Some repNames correspond to multiple repeat familys
    #' Split those repNames with multiple repeat families such that
    #' each repName correspond to a single repeat family.
    #' those include (CATTC)n, (GAATC)n, (GAATG)n, MamRep137,
    #' MamRep1879, MER132, MER96, MER96B
    if (finerSplitByRepFamily)
        features <- .finerSplitByRepeatFamily(features)

    #' get element's repClass and repFamily to metadata
    meta_data <- .getElementMetaData(features)
    mcols(features) <- meta_data
    return(features)
}
