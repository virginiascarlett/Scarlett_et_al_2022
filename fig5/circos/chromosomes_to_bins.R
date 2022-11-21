# The first script for producing a histogram track in circos that will 
# give each bin either a 0 or a 1, useful if you just want to indicate
# the presence/absence of a feature. Follow this script with make_01_histo_circos.R,
# but you will also need to run prepare_centro_files_circos.R first.
# Run with genome name as an argument on the command line like so:
# Rscript chromosomes_to_bins.R ABR113

options(scipen=999) #to prevent scientific notation

args = commandArgs(trailingOnly=TRUE)
genome <- args[1]

setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')

#chrdat <- read.table('karyotype.Bd21Osat.txt', stringsAsFactors=FALSE)
chrdat <- read.table(sprintf('karyotype.%s.txt', genome), stringsAsFactors=FALSE)
chrdat <- chrdat[,c(4,6)]

get_chrom_df <- function(myrow){
  starts <- seq(1,chrdat[myrow,2],by=100000)
  ends <- append(sapply(starts, function(x) x-1)[-1], chrdat[myrow,2])
  chromdf <- data.frame(rep(chrdat[myrow,1], times=length(starts)), starts, ends)
  return(chromdf)
}

finaldf <- do.call("rbind", lapply(seq(nrow(chrdat)), function(x) get_chrom_df(x)))
write.table(finaldf, file=sprintf("bins.%s.txt", genome), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=",")
