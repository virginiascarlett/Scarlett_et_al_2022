library(GenomicRanges)
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
  strand=frags$strand)


m <- mergeByOverlaps(fragGR, binGR)
temp <- m[,2]
temp$TEspace <- width(m[,1]) #a df where each row represents a TE fragment, and the length of the TE fragment is a metadata column
q <- findOverlaps(binGR, temp, type='equal') #"Specifying equal as the type returns the intersection of the start and end matches"
#queryHits of q is a genomic bin (row of binGR), subjectHits is a TE fragment (row of temp)
indices <- queryHits(q)
binGR$TEspace <- rep(0, times=length(binGR))
#takes about a minute
j <- sapply(unique(indices), function(i) binGR[i]$TEspace <<- sum(temp[indices==i]$TEspace))
final <- data.frame(Chr=seqnames(binGR), Start=start(binGR), End=end(binGR), TEDensity=round(binGR$TEspace/width(binGR), 3))

write.table(final, file=sprintf("TEdensity.%s.%dkb.tsv", genome, windowSize/1000), row.names = FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(final, file=sprintf("TEdensity.%s.%dkb.csv", genome, windowSize/1000), row.names = FALSE, col.names=TRUE, quote=FALSE, sep=",")















