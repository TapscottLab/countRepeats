# Overview

The _countRepeats_ package is a R/Bioconductor package providing tools to count hits for repeat elements. The method is adapted from GenomicAlignments::summarizeOverlaps() with a few adjustments. A read must overlap at least 20 bps with the repeat elements to be counted, and the count is adjusted by dividing the number of reported multiple alignments (NH). 

Using tools from the DESeq2 and ggplot2 packages, _countRepeats_ provides downstream analysis including enrichment/depletion analysis and visualization.
