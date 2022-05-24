suppressPackageStartupMessages({
  library(GenomicRanges)
})
options(scipen=999) #to prevent scientific notation

get_chrom_df <- function(myrow){
  starts <- seq(1,chrdat[myrow,2],by=windowSize)
  ends <- append(sapply(starts, function(x) x-1)[-1], chrdat[myrow,2])
  chromdf <- data.frame(rep(chrdat[myrow,1], times=length(starts)), starts, ends)
  return(chromdf)
}

############# EDIT HERE ################
setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')
genome <- 'Bhyb26'
chromsizesfile <- paste0(genome, '.chrom.sizes')
fragfile <- paste0('allfragments.classified_', genome)
windowSize <- 200000
########################################


chrdat <- read.csv(chromsizesfile, sep='\t', header=FALSE, stringsAsFactors = FALSE)
frags <- read.csv(fragfile, header=TRUE, sep='\t', stringsAsFactors = FALSE)

#A GRanges object of the whole genome divided into <windowSize>-sized bins
bindf <- do.call("rbind", lapply(seq(nrow(chrdat)), function(x) get_chrom_df(x)))
binGR <- GRanges(
  seqnames = bindf[,1],
  ranges = IRanges(bindf[,2], bindf[,3]),
)

fragGR <- GRanges(
  seqnames=frags$chrom, 
  ranges=IRanges(start=frags$frag_start, end=frags$frag_end),
  strand=frags$strand, 
  exemplar=frags$exemplar,
  classification=frags$classification)

ov <- findOverlaps(fragGR, binGR)
#takes like 30 seconds
binGR$numfams <- sapply(seq(length(binGR)), function(i) length(unique(fragGR[queryHits(subset(ov, subjectHits==i)), ]$exemplar)) )
#takes like 30 seconds
binGR$numRTfams <- sapply(seq(length(binGR)), function(i) length(unique(
  subset(fragGR[queryHits(subset(ov, subjectHits==i)), ], grepl( "retrotransposon", classification, fixed = TRUE))$exemplar)) )

binGR$numDNAfams <- sapply(seq(length(binGR)), function(i) length(unique(
  subset(fragGR[queryHits(subset(ov, subjectHits==i)), ], grepl( "DNA_transposon", classification, fixed = TRUE))$exemplar)) )

final <- as.data.frame(binGR)[,c(1, 2, 3, 6, 7, 8)]
names(final) <- c("Chromosome", "Start", "End", "Num_Fams", "Num_RT_Fams", "Num_DNA_Fams")
write.table(final, file=sprintf("TEdiversity.%s.%dkb.csv", genome, windowSize/1000), row.names = FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(final[,c(1:4)], file=sprintf("TEdiversity.%s.%dkb.circos.tsv", genome, windowSize/1000), row.names = FALSE, col.names=FALSE, quote=FALSE, sep="\t")

