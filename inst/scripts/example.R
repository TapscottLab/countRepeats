#'test
library(GenomicAlignments)
library(canFam3.rmsk)
library(rmskStats)
canFam3.rmsk <- keepStandardChromosomes(canFam3.rmsk, pruning.mode="coarse")
features <- splitByRepName(canFam3.rmsk)

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "all.reads.rda"))
se <- rmskStats::summarizeOverlaps.adjNH(features=features,
                                              all.reads=all.reads[[1]],
                                              type="within",
                                              ignore.strand=TRUE,
                                              inter.feature=TRUE)
rowData(se)

## give an example for counting HSATII
