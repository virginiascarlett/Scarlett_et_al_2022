suppressPackageStartupMessages({
  library(GenomicRanges)
  library(stringr)
  library(ggplot2)
})
options(scipen=999) #to prevent scientific notation

#This script checks overlap between genes and TEs. Input is a gff that combines the phytozome gene annotation 
#and my own TE annotation with the INDIVIDUAL TE library. The latter was produced by my script prep_McClintock_raw_gff.py.
#Then the two gffs were concatenated like so:
#awk 'NR > 3 { print }' Bhybridum_463_v1.1.gene.gff3 > ABR113.noheader.gff 
#cat ABR113_RepeatMasker_TElib.gff ABR113.noheader.gff > ABR113_gff_w_repeats.new.gff
#bedtools sort -i ABR113_gff_w_repeats.new.gff > ABR113_gff_w_repeats.new.sorted.gff

#awk 'NR > 6 { print }' BhybridumBhyb26v2.1.gene.gff3 > Bhyb26.noheader.gff 
#cat Bhyb26_RepeatMasker_TElib.gff Bhyb26.noheader.gff > Bhyb26_gff_w_repeats.new.gff
#bedtools sort -i Bhyb26_gff_w_repeats.new.gff > Bhyb26_gff_w_repeats.new.sorted.gff

#The middle portion requires a file with coordinates for centromeric and pericentromeric regions,
#i.e. Bhyb26centromeres.txt as well as a file with the size of each chromosome i.e. Bhyb26.chrom.sizes.

#The last portion requires a 'pseudogeneIntervals.txt' file for the lost genes.

gff2GR <- function(mygffdf) {
  return(GRanges(
    seqnames=mygffdf[,1], 
    ranges=IRanges(start=mygffdf[,4], end=mygffdf[,5]),
    strand=mygffdf[,7],
    feature=mygffdf[,3],
    annotsource=mygffdf[,2],
    info=mygffdf[,9])
  )
}
get_centro_GR <- function(genome) {
  centrofile <- paste0(genome, 'centromeres.txt') #produced manually and painstakingly
  # Creating unequal genomic bins for centromeric, pericentromeric, and "distal" regions
  centrodat <- read.csv(centrofile, sep='\t', header=FALSE, stringsAsFactors = FALSE)
  names(centrodat) <- c('chromosome', 'start', 'end', 'feature')
  centrodat$start <- centrodat$start*1000000
  centrodat$end <- centrodat$end*1000000
  centroGR <- GRanges(  
    seqnames = centrodat$chromosome,
    ranges = IRanges(centrodat$start, centrodat$end),
    feature = centrodat$feature
  )
  #In my input file, pericentromeres encompass the centromeres. Let's get them as non-overlapping features now.
  non_overlapping <- disjoin(centroGR)
  just_centros <- subset(centroGR, feature=='centromere')
  non_overlapping$feature <- rep('*', times=length(non_overlapping))
  q <- findOverlaps(non_overlapping, just_centros, type='equal')
  non_overlapping[queryHits(q),]$feature <- 'centromere'
  non_overlapping[which(non_overlapping$feature=='*'),]$feature <- 'pericentromere'
  centroGR <- non_overlapping #just renaming it
  return(centroGR)
}
get_distal <- function(cenGR, mychr, chrdat){
  tempGR <- GRanges(
    seqnames=mychr,
    ranges=IRanges(1, start(cenGR[min(which(seqnames(cenGR)==mychr)),])-1),
    strand='*',
    feature='distal')
  tempGR2 <- GRanges(
    seqnames=mychr,
    ranges=IRanges(end(cenGR[max(which(seqnames(cenGR)==mychr)),])+1, chrdat[which(chrdat$V1==mychr),"V2"]),
    strand='*',
    feature='distal')
  return(c(tempGR,tempGR2))
}
get_overlap_within_feature <- function(geneGR, TEGR, featureGR, myfeature){
  genes_sub <- geneGR[queryHits(findOverlaps(geneGR, subset(featureGR, feature==myfeature))),]
  TEs_sub <- TEGR[queryHits(findOverlaps(TEGR, subset(featureGR, feature==myfeature))),]
  ov_sub <- findOverlaps(genes_sub, TEs_sub, ignore.strand=F) #NOTE STRAND CHOICE HERE
  return(length(unique(queryHits(ov_sub)))/length(genes_sub)) #what percent of genes overlap a TE IN e.g. DISTAL REGIONS?
}
#setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')



########## PART 1: GET OVERLAP BETWEEN TEs AND GENES ########## 
#note to self: maybe check LTR-RTs separately

Bhyb26 <- read.csv('/home/virginia/Documents/school/vogelLab/Bhyb26_gff_w_repeats.new.sorted.gff', sep="\t", header=FALSE, stringsAsFactors = FALSE)
ABR113 <- read.csv('/home/virginia/Documents/school/vogelLab/ABR113_gff_w_repeats.new.sorted.gff', sep="\t", header=FALSE, stringsAsFactors = FALSE)
ABR114 <- read.csv('/home/virginia/Documents/school/vogelLab/ABR114_gff_w_repeats.new.sorted.gff', sep="\t", header=FALSE, stringsAsFactors = FALSE)
Bd21 <- read.csv('/home/virginia/Documents/school/vogelLab/Bd21_gff_w_repeats.new.sorted.gff', sep="\t", header=FALSE, stringsAsFactors = FALSE)

Bhyb26GR <- gff2GR(Bhyb26)
ABR113GR <- gff2GR(ABR113)
ABR114GR <- gff2GR(ABR114)
Bd21GR <- gff2GR(Bd21)
Bhyb26D <- Bhyb26GR[which(startsWith(as.character(seqnames(Bhyb26GR)), "BhD"))]
Bhyb26S <- Bhyb26GR[which(startsWith(as.character(seqnames(Bhyb26GR)), "BhS"))]
ABR113D <- ABR113GR[which(startsWith(as.character(seqnames(ABR113GR)), "BhD"))]
ABR113S <- ABR113GR[which(startsWith(as.character(seqnames(ABR113GR)), "BhS"))]

gene_TE_overlap <- function(myGR) {
  mygenes <- subset(myGR, feature=="gene")
  myTEs <- subset(myGR, annotsource =="RepeatMasker")
  myov <- findOverlaps(mygenes, myTEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
  return( length(unique(queryHits(myov)))/length(mygenes)*100 )#what percent of genes overlap a TE?
}
exon_TE_overlap <- function(myGR) {
  myCDS <- subset(myGR, feature=="CDS")
  myTEs <- subset(myGR, annotsource =="RepeatMasker")
  myov <- findOverlaps(myCDS, myTEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
  return( length(unique(queryHits(myov)))/length(myCDS)*100 )#what percent of exons overlap a TE?
}
LTRRT_overlap <- function(myGR) {
  mygenes <- subset(myGR, feature=="gene")
  myTEs <- subset(myGR, feature == "LTR_retrotransposon_RLG" | feature == "LTR_retrotransposon_RLC" | feature == "LTR_retrotransposon_RLX")
  myov <- findOverlaps(mygenes, myTEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
  return( length(unique(queryHits(myov)))/length(mygenes)*100 )#what percent of genes overlap a TE?
}
distance2TEs <- function(myGR) {
  mygenes <- subset(myGR, feature=="gene")
  myTEs <- subset(myGR, annotsource =="RepeatMasker")
  myov <- findOverlaps(mygenes, myTEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
  noov <- mygenes[-unique(queryHits(myov))]
  dhits <- distanceToNearest(noov, myTEs, ignore.strand=F)
  mydistances <- as.data.frame(dhits)$distance
  return( mean(mydistances) )
}

#Percent of Bhyb26 genes that overlap a TE:
gene_TE_overlap(Bhyb26D)
#Percent of Bhyb26 exons that overlap a TE:
exon_TE_overlap(Bhyb26D)
#Percent of Bhyb26 genes that overlap a LTR-RT:
LTRRT_overlap(Bhyb26D)
#Distance to nearest TE, Bhyb26 genes:
distance2TEs(Bhyb26D)

#Percent of Bhyb26 genes that overlap a TE:
gene_TE_overlap(Bhyb26S)
#Percent of Bhyb26 exons that overlap a TE:
exon_TE_overlap(Bhyb26S)
#Percent of Bhyb26 genes that overlap a LTR-RT:
LTRRT_overlap(Bhyb26S)
#Distance to nearest TE, Bhyb26 genes:
distance2TEs(Bhyb26S)

#Percent of ABR113 genes that overlap a TE:
gene_TE_overlap(ABR113D)
#Percent of ABR113 exons that overlap a TE:
exon_TE_overlap(ABR113D)
#Percent of ABR113 genes that overlap a LTR-RT:
LTRRT_overlap(ABR113D)
#Distance to nearest TE, ABR113 genes:
distance2TEs(ABR113D)

#Percent of ABR113 genes that overlap a TE:
gene_TE_overlap(ABR113S)
#Percent of ABR113 exons that overlap a TE:
exon_TE_overlap(ABR113S)
#Percent of ABR113 genes that overlap a LTR-RT:
LTRRT_overlap(ABR113S)
#Distance to nearest TE, ABR113 genes:
distance2TEs(ABR113S)

#Percent of ABR114 genes that overlap a TE:
gene_TE_overlap(ABR114GR)
#Percent of ABR114 exons that overlap a TE:
exon_TE_overlap(ABR114GR)
#Percent of ABR114 genes that overlap a LTR-RT:
LTRRT_overlap(ABR114GR)
#Distance to nearest TE, ABR114 genes:
distance2TEs(ABR114GR)

#Percent of Bd21 genes that overlap a TE:
gene_TE_overlap(Bd21GR)
#Percent of Bd21 exons that overlap a TE:
exon_TE_overlap(Bd21GR)
#Percent of Bd21 genes that overlap a LTR-RT:
LTRRT_overlap(Bd21GR)
#Distance to nearest TE, Bd21 genes:
distance2TEs(Bd21GR)



#Let's just keep the genes and the TEs (intronic and UTR TEs still count)
Bhyb26genes <- subset(Bhyb26GR, feature=="gene")
Bhyb26TEs <- subset(Bhyb26GR, annotsource =="RepeatMasker")
Bhyb26ov <- findOverlaps(Bhyb26genes, Bhyb26TEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
print("Percent of Bhyb26 genes that overlap a TE:")
print( length(unique(queryHits(Bhyb26ov)))/length(Bhyb26genes) )#what percent of genes overlap a TE?

ABR113genes <- subset(ABR113GR, feature=="gene")
ABR113TEs <- subset(ABR113GR, annotsource =="RepeatMasker")
ABR113ov <- findOverlaps(ABR113genes, ABR113TEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
print("Percent of ABR113 genes that overlap a TE:")
print( length(unique(queryHits(ABR113ov)))/length(ABR113genes) )

#WHAT ABOUT JUST EXONS?
Bhyb26CDS <- subset(Bhyb26GR, feature=="CDS")
Bhyb26ovtmp <- findOverlaps(Bhyb26CDS, Bhyb26TEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
print( length(unique(queryHits(Bhyb26ovtmp)))/length(Bhyb26CDS) )#what percent of exons overlap a TE?

ABR113CDS <- subset(ABR113GR, feature=="CDS")
ABR113ovtmp <- findOverlaps(ABR113CDS, ABR113TEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
print( length(unique(queryHits(ABR113ovtmp)))/length(ABR113CDS) )#what percent of exons overlap a TE?





#DISTANCE TO NEAREST TE (FOR THOSE GENES THAT DON'T OVERLAP A TE)?
#get genes that don't overlap a TE
Bhyb26noov <- Bhyb26genes[-unique(queryHits(Bhyb26ov))]
Bhyb26dhits <- distanceToNearest(Bhyb26noov, Bhyb26TEs, ignore.strand=F)
Bhyb26distances <- as.data.frame(Bhyb26dhits)$distance
mean(Bhyb26distances)

ABR113noov <- ABR113genes[-unique(queryHits(ABR113ov))]
ABR113dhits <- distanceToNearest(ABR113noov, ABR113TEs, ignore.strand=F)
ABR113distances <- as.data.frame(ABR113dhits)$distance
mean(ABR113distances)
###################################################################

######### PART 2: CHECK GENE/TE OVERLAP FOR CENTROMERIC, PERICENTROMERIC, AND DISTAL GENES ######### 

Bhyb26cen <- get_centro_GR('Bhyb26')
Bhyb26chrdat <- read.csv('Bhyb26.chrom.sizes', sep='\t', header=FALSE, stringsAsFactors = FALSE)
seqlevels(Bhyb26cen) <- Bhyb26chrdat$V1
Bhyb26cen <- sort(Bhyb26cen) 
L_of_GRs <- sapply(Bhyb26chrdat$V1, function(mychr) get_distal(Bhyb26cen, mychr, Bhyb26chrdat))
Bhyb26distalGR <- unlist(as(L_of_GRs, "GRangesList"))
Bhyb26features <- c(Bhyb26cen, Bhyb26distalGR)
seqlevels(Bhyb26features) <- Bhyb26chrdat$V1
Bhyb26features <- sort(Bhyb26features) 

ABR113cen <- get_centro_GR('ABR113')
ABR113chrdat <- read.csv('ABR113.chrom.sizes', sep='\t', header=FALSE, stringsAsFactors = FALSE)
seqlevels(ABR113cen) <- ABR113chrdat$V1
ABR113cen <- sort(ABR113cen) 
L_of_GRs <- sapply(ABR113chrdat$V1, function(mychr) get_distal(ABR113cen, mychr, Bhyb26chrdat))
ABR113distalGR <- unlist(as(L_of_GRs, "GRangesList"))
ABR113features <- c(ABR113cen, ABR113distalGR)
seqlevels(ABR113features) <- ABR113chrdat$V1
ABR113features <- sort(ABR113features) 

print("Percent of Bhyb26 genes that overlap a TE, distal, pericentromere, and centromere regions only:")
print( get_overlap_within_feature(Bhyb26genes, Bhyb26TEs, Bhyb26features, "distal") )
print( get_overlap_within_feature(Bhyb26genes, Bhyb26TEs, Bhyb26features, "pericentromere") )
print( get_overlap_within_feature(Bhyb26genes, Bhyb26TEs, Bhyb26features, "centromere") )
print("Percent of ABR113 genes that overlap a TE, distal, pericentromere, and centromere regions only:")
print( get_overlap_within_feature(ABR113genes, ABR113TEs, ABR113features, "distal") )
print( get_overlap_within_feature(ABR113genes, ABR113TEs, ABR113features, "pericentromere") )
print( get_overlap_within_feature(ABR113genes, ABR113TEs, ABR113features, "centromere") )
############################################################## 



############################### PART 3: CHECK GENE/TE OVERLAP FOR MY 'LOST GENES' ############################### 
#WHAT ABOUT 
#Now that I am no longer using parseAlignmentsFromExonerateNEWNEW.py as part of my pipeline, 
#that script had to be fixed up. Ran on my 'lost genes'. For the conserved genes,
#just get the start and end of gene from gff.

lost <- read.csv('/home/virginia/Documents/school/vogelLab/notebook/2021/lost_genes/pseudogeneIntervals.dipsvsBhyb26.txt', sep='\t', header=F, stringsAsFactors = F)
lostGR <- GRanges(
  seqnames=lost[,"V6"],
  ranges=IRanges(start=lost[,"V7"], end=lost[,"V8"]),
  strand=lost[,"V9"],
  dipgene=lost[,"V1"])
lostGR_condensed <- unlist(range(split(lostGR, ~dipgene))) #just take the start and end of the alignable region
#if there is a TE insertion, we don't expect that to be alignable
Bhyb26LGov <- findOverlaps(lostGR_condensed, Bhyb26TEs, ignore.strand=F) #NOTE STRAND CHOICE HERE
length(unique(queryHits(Bhyb26LGov))) #this is the number of lost genes that overlap a TE (130)

#Get the number of genes that overlap a TE, for all trials

conservedgenes <- read.csv('/home/virginia/Documents/school/vogelLab/notebook/2021/lost_genes/conservedgenes.txt', sep='\t', header=T, stringsAsFactors = F) #produced by my script bin_orthogroups.R
conservedgenes <- c(conservedgenes$Bhyb26Dgene, conservedgenes$Bhyb26Sgene)
Bhyb26genes$geneID <- gsub("ID=", "", sapply(strsplit(Bhyb26genes$info, ";Name="), "[[", 1))
Bhyb26genesslim <- subset(Bhyb26genes, geneID %in% conservedgenes)
trialinfo <- read.csv('/home/virginia/Documents/school/vogelLab/notebook/2021/lost_genes/picks_by_trial.txt', sep='\t', header=T, stringsAsFactors = F) #produced by my script gather_slurmIDs.py

Bhyb26ovGR <- Bhyb26genes[unique(queryHits(Bhyb26ov))]
Bhyb26ov_genes <- gsub("ID=", "", sapply(strsplit(Bhyb26ovGR$info, ";Name="), "[[", 1))
trialinfo$TEov <- trialinfo$Hybgene %in% Bhyb26ov_genes
res <- rbind(by(trialinfo$TEov, trialinfo$Trial, sum))
mean(res)
# [1] 209.585
#^On average, my control trials with conserved genes yielded 209.85 genes that overlap a TE





#'Vox referee vox dei': reviewer requested the distribution of outcomes (Aug 2022)
toPlot<-as.data.frame(t(data.frame('outcome'=res)))
names(toPlot) <- c('outcome')
p <- ggplot(toPlot, aes(x=outcome))+
  geom_histogram(color = "black", fill="#56B4E9") + 
  geom_vline(xintercept=130, color='red', size=0.7, linetype='dashed')+
  ylab("Number of trials")+
  xlab("Number of genes overlapping a TE")+ 
  theme(
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/aug12_TEoverlap_trials.tiff", units="in", width=5, height=3, res=300)
p+ 
  #scale_y_continuous(expand = c(0, 0), limits = c(0, max(ggplot_build(p3)$data[[1]]$ymax)*1.1))+
  scale_x_continuous(expand = c(0, 0), limits = c(110, 260))
dev.off()







