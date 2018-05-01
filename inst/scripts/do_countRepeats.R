#!/app/easybuild/software/R/3.4.3-foss-2016b-fh1/bin/Rscript
#./do_countHSATII.R
#SBATCH -M beagle -p largenode --exclusive -n1 

#' This R script is an example for counting HSATII on hg38 genome. 
#' supports embarrassingly parallel programming (called by example_epp.R)
#' by taking one sample at a time and run on largenode cluster.
#' Note that the parameter of countRepeats() is set specifically for
#' a particular bam files with paired-end reads in which the strand of
#' origin is carried out by the second reads. One should make the
#'  parameters to fit his/her circumstences.
#' 
#' Input:
#' -b, --bamfile  bam file name
#' -r, --repeats  repeat features file name
#' -o, --out      output directory
#'
#' Get help:
#' Rscript --vanilla do_countHSATII.R --help

#' parse argument using python style
library(optparse)
option_list <- list(make_option(c("-b", "--bamfile"), type="character",
                                default=NULL,
                                help="bam file name", metavar="character"),
                    make_option(c("-r", "--repeats"), type="character",
                                default=NULL,
                                help="repeat features file name",
                                metavar="character"),
                    make_option(c("-o", "--out"), type="character",
                                default=NULL,
                                help="output directory",
                                metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#' sanity check
if (!file.exists(opt$repeats)) stop("Repeat file does not exist")
if (!file.exists(opt$bamfile)) stop("The bam file does not exist")

#' count start - strand aware, second reads carry strand of origin
library(countRepeats)
library(GenomicAlignments)
features <- get(load(opt$repeats))
features <- keepStandardChromosomes(features, pruning.mode="coarse",
                                    species="Homo_sapiens")

sample_name <- sub(".bam", "", basename(opt$bamfile))
se <- countRepeats(features=features,
                   bam_files= opt$bamfile,
                   singleEnd=FALSE,
                   type="any",
                   ignore.strand=FALSE,
                   inter.feature=FALSE,
                   round.up=FALSE,
                   strandMode=2)
colnames(se) <- sample_name
rownames(se) <- features$exon_name
save(se, file=file.path(opt$out, paste0(sample_name, "_hsat2.rda")))
     
