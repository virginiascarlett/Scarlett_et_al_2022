# This is a newer version of make_circos_KEE_histo.R.
# Run with bin file (produced by chromosomes_to_bins.R) and your regions 
# of interest (3 columns, space-delimited) produced by prepare_centro_file_circos.R
# as arguments on the command line like so:
# Rscript make_01histo_circos.R bins.ABR113.txt ABR113centromeres.circos.txt
# Do for both pericentromeres and centromeres!

suppressPackageStartupMessages({
  library(GenomicRanges)
})

setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')

args = commandArgs(trailingOnly=TRUE)
binfile <- args[1]
featfile <- args[2]
  
bindf <- read.csv(binfile, stringsAsFactors = FALSE, header=FALSE)
featdf <- read.csv(featfile, stringsAsFactors = FALSE, header=FALSE, sep=" ")

binGR <- GRanges(
  seqnames = bindf[,1],
  ranges = IRanges(bindf[,2], bindf[,3]),
  score = 0
)
featGR <- GRanges(
  seqnames = featdf[,1],
  ranges = IRanges(featdf[,2], featdf[,3])
)
m <- findOverlaps(binGR, featGR)
binGR[queryHits(m),]$score <- 1

final <- as.data.frame(binGR)[,c(1,2,3,6)]

write.table(final, file=sprintf("%s.circos.final.txt", strsplit(featfile, "\\.")[[1]][1]), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")


