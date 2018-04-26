#' ml R/3.4.0-foss-2016b-fh1
#' try hg19
library(BSgenome.Hsapiens.UCSC.hg19)

library(GenomicAlignments)
library(Biostrings)
library(rmskStats)
library(ShortRead)
library(parallel)

#' input
BSgenome <- BSgenome.Hsapiens.UCSC.hg19
rmsk <- get(load("~/tapscott/hg19/hg19.rmsk.rda"))
outDir <- file.path("~/tapscott/hg19")
species <- "Homo_sapiens" #names(genomeStyles())
cores <- 4L
keep.standard.chr <- TRUE
name <- "hg19.RMSK.repName.fa"

#' (1) make sure the load regions are  not beyond the boundaries of non-circular sequence
name <- "hg19.RMSK.repName.fa"
fa_filename <- file.path(outDir, name)
seq_info <- seqinfo(BSgenome)
seq_length <- seqlengths(seq_info)
rmsk <- keepStandardChromosomes(rmsk, species=species,
                                pruning.mode="coarse")
rmsk <- rmskStats::splitByRepName(rmsk=rmsk)

#' 
#' (2) get sequence for each RepName
#' For each repName, concatenate the sequence of the entities and separated by "N-150" customized intron
#' ("NNN-NNNN" with 150 width)
#' 
.detectFalseStartSite <- function(gr) {
    #' is there any start site at 0 -> correct to 1
    idx <- start(gr) == 0
    if (sum(idx) > 0) start(gr)[idx] = 1
    gr
}

.detectFasleEndSite <- function(gr) {
}


intron <- strrep("N", 150)
repName_seq <- mclapply(names(rmsk), function(x) {
    print(x)
    gr <- rmsk[[x]]
    gr <- .detectFalseStartSite(gr)
    #'gr <- .detectFalseEndSite(gr, BSgenome)
    
    seq <- getSeq(BSgenome, gr)
    tmp <- as.character(seq)
    tmp <- paste0(tmp, collapse=intron)
    seq <- DNAString(paste0(intron, tmp, intron))
}, mc.cores=cores, mc.preschedule=FALSE)
names(repName_seq) <- names(rmsk)
fa <- DNASTringSet(unlist(repName_seq))
writeFasta(fa, file=fa_filename)


## debuging
gr <- rmsk[["L1_Rat1"]]
a <- gr[seqnames(gr)=="chrX"]

for (i in 1:length(a)) {
    print(i)
    seq <- getSeq(BSgenome.Rnorvegicus.UCSC.rn6, a[i])
    ## conclusion: a[20] has the problem, the start is 0
    ##GRanges object with 1 range and 8 metadata columns:
    ##     seqnames    ranges strand |   genoLeft  repName repClass repFamily
    ##        <Rle> <IRanges>  <Rle> |  <integer> <factor> <factor>  <factor>
    ##[1]      chrX  [0, 668]      + | -159969353  L1_Rat1     LINE        L1
    ##        repStart    repEnd   repLeft        id
    ##        <integer> <integer> <integer> <integer>
    ##[1]      5472      6142      -486         4
}
