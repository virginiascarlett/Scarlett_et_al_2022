# Prepares non-overlapping centromere and pericentromere files for 
# making a histogram track in circos. Runthis after chromosomes_to_bins.R
# and follow this script with make_01histo_circos.R.
# Run with genome name as an argument on the command line like so:
# Rscript prepare_centro_files_circos.R ABR113
suppressPackageStartupMessages({
  library(GenomicRanges)
})

options(scipen=999) #to prevent scientific notation
get_centro_data <- function(filename) {
  dat <- read.csv(filename, stringsAsFactors = FALSE, header=FALSE, sep="\t")
  dat[,c(2, 3)] <- dat[,c(2,3)]*1000000
  cen <- subset(dat, V4=="centromere")
  peri <- subset(dat, V4=="pericentromere")
  return(list(cen, peri))
}

args = commandArgs(trailingOnly=TRUE)
genome <- args[1]

setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')

datL <- get_centro_data(sprintf("%scentromeres.txt", genome))
cendat <- datL[[1]]
perdat <- datL[[2]]
cenGR <- GRanges(
  seqnames = cendat[,1],
  ranges = IRanges(cendat[,2], cendat[,3])
)
perGR <- GRanges(
  seqnames = perdat[,1],
  ranges = IRanges(perdat[,2], perdat[,3])
)
perfinal <- setdiff(perGR, cenGR)

write.table(as.data.frame(cenGR)[,c(1,2,3)], file=sprintf("%scentromeres.circos.txt", genome), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")
write.table(as.data.frame(perfinal)[,c(1,2,3)], file=sprintf("%spericentromeres.circos.txt", genome), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")

