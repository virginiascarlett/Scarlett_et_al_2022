#### Get subgenome-specific TEs as a "density" file for my script frag_heatmap.R ###
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(scales)
})
options(scipen=999) #to prevent scientific notation

get_chrom_df <- function(myrow){
  starts <- seq(1,chrdat[myrow,2],by=windowSize)
  ends <- append(sapply(starts, function(x) x-1)[-1], chrdat[myrow,2])
  chromdf <- data.frame(rep(chrdat[myrow,1], times=length(starts)), starts, ends)
  return(chromdf)
}

get_specificity <- function(numD, numS){
  return(round(
    max( c(numD/sum(numD, numS), numS/sum(numD, numS))
    ), 3))
}
# ^ Gets the proportion of TEs on each subgenome and returns the larger number.
# For example, Copia.chain00774#LTR/Copia          788         1122
# Returns 0.587. That means that 41% of TE copies are on one
# subgenome and 59% are on the other subgenome.


setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')
genome <- 'ABR113'
chromsizesfile <- paste0(genome, '.chrom.sizes')
fragfile <- paste0('allfragments.classified_', genome) #, '_masterannot')
subgspecfile <- paste0('subg_specificity_by_family.txt_', genome) #, '_master') #from step7 of my TE pipeline
windowSize <- 200000
fragdat <- read.csv(fragfile, sep='\t', header=TRUE, stringsAsFactors = FALSE)
fragGR <- GRanges(
  seqnames=fragdat[,2], 
  ranges=IRanges(start=fragdat[,3], end=fragdat[,4]),
  strand=fragdat[,9],
  TEcode = fragdat[,6],
  exemplar = fragdat[,10]
)
subg_all <- read.csv(subgspecfile, sep='\t', header=TRUE, stringsAsFactors = FALSE)
subg_all$specificity <- sapply(seq(nrow(subg_all)), function(n) get_specificity(subg_all[n,2], subg_all[n,3]) )

# Defining subgenome-specific TE families as those with more than 90% of copies on one subgenome, and 5 or more copies total
#ADDED NEW CODE AUG 25 2021, I HAVE NO IDEA WHY THIS DOESN'T WORK:
#subgspecfams <- subset(subg_all, specificity >= 0.90 & sum(num_copies_D, num_copies_S) >= 5)
#IT GETS THE SPECIFICITY PART BUT IT IGNORES THE FAMILY SIZE PART
subgspecfams <- subset(subg_all, specificity >= 0.90)
subgspecfams <- subgspecfams[mapply(function(x, y) x+y >= 5, subgspecfams$num_copies_D, subgspecfams$num_copies_S),]
subgspecGR <- subset(fragGR, exemplar %in% subgspecfams$famname)

chrdat <- read.csv(chromsizesfile, sep='\t', header=FALSE, stringsAsFactors = FALSE)
bindf <- do.call("rbind", lapply(seq(nrow(chrdat)), function(x) get_chrom_df(x)))
binGR <- GRanges(
  seqnames = bindf[,1],
  ranges = IRanges(bindf[,2], bindf[,3]),
)
binGR$score <- rep(0, times=length(binGR))

m <- findOverlaps(binGR, subgspecGR)
binGR[as.numeric(names(table(queryHits(m)))),]$score <- unlist(unname(table(queryHits(m))))
final <- as.data.frame(binGR)
final <- final[,c(1,2,3,6)]
names(final) <- c("Chr", "Start", "End", "SubgSpecTEs")
write.table(final, file=paste0('subgenome_specific_TEs.', genome, '.csv'), sep=",", quote=FALSE, row.names=FALSE)
write.table(final, file=paste0('subgenome_specific_TEs.', genome, '.circos.tsv'), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
#write.table(final, file=paste0('subgenome_specific_TEs.Bhyb26.80perc.csv'), sep=",", quote=FALSE, row.names=FALSE)

################################################################3
#Adding some stuff, Aug 10 2021

# favorsD <- subset(subg_all, num_copies_D > num_copies_S)
# favorsD$domsubg <- 'D'
# favorsD$famsize <- mapply(function(x, y) sum(x, y), favorsD$num_copies_D, favorsD$num_copies_S)
# favorsS <- subset(subg_all, num_copies_S > num_copies_D)
# favorsS$domsubg <- 'S'
# favorsS$famsize <- mapply(function(x, y) sum(x, y), favorsS$num_copies_D, favorsS$num_copies_S)
# 
# toPlot <- rbind(favorsD, favorsS)
# toPlot <- subset(toPlot, famsize>=5)

#dim(subset(toPlot, domsubg=='D' & specificity<0.6))
#dim(subset(toPlot, domsubg=='S' & specificity<0.6))


ABR113 <- final
ABR113_fragGR <- fragGR
TEgr_ABR113 <- subset(ABR113_fragGR, !duplicated(TEcode))
TEgr_ABR113$isspecific <- TEgr_ABR113$exemplar %in% subgspecfams$famname
#CRITICALLY IMPORTANT:
#now go back and re-run the above code with the other genome!!!!
#Run this chunk of code^ or the one below depending on your current genome
Bhyb26 <- final
Bhyb26_fragGR <- fragGR
TEgr_Bhyb26 <- subset(Bhyb26_fragGR, !duplicated(TEcode))
TEgr_Bhyb26$isspecific <- TEgr_Bhyb26$exemplar %in% subgspecfams$famname

#now we can proceed.
#% of TEs that come from subgenome-specific families
sum(TEgr_ABR113$isspecific)/length(TEgr_ABR113)
sum(TEgr_Bhyb26$isspecific)/length(TEgr_Bhyb26)
    

ABR113_by_chr <- sapply(unique(ABR113$Chr), function(mychr) sum(ABR113[which(ABR113$Chr==mychr),]$SubgSpecTEs))
Bhyb26_by_chr <- sapply(unique(Bhyb26$Chr), function(mychr) sum(Bhyb26[which(Bhyb26$Chr==mychr),]$SubgSpecTEs))
bychr_long <- data.frame(genome=c(rep("ABR113", times=15), rep("Bhyb26", times=15)), 
                         chrom=c(unique(as.character(ABR113$Chr)), unique(as.character(Bhyb26$Chr))), 
                         num_subg_spec_TEs=c(ABR113_by_chr, Bhyb26_by_chr), 
                         stringsAsFactors = F)

totals <- rbind( data.frame(
  genome = "ABR113",
  chrom = unique(seqnames(TEgr_ABR113)),
  total_TEs = sapply(unique(seqnames(TEgr_ABR113)), function(mychr) length(subset(TEgr_ABR113, seqnames==mychr))),
  stringsAsFactors = FALSE),
  data.frame(
    genome = "Bhyb26",
    chrom = unique(seqnames(TEgr_Bhyb26)),
    total_TEs = sapply(unique(seqnames(TEgr_Bhyb26)), function(mychr) length(subset(TEgr_Bhyb26, seqnames==mychr))),
    stringsAsFactors = FALSE)
)

bychr_totals <- merge(bychr_long, totals, by.x=c('genome', 'chrom'), by.y=c('genome','chrom'))
bychr_totals$perc <- bychr_totals$num_subg_spec_TEs/bychr_totals$total_TEs
bychr_totals$chrom <- factor(bychr_totals$chrom, levels = c("BhD1", "BhD2",  "BhD3",  "BhD4",  "BhD5",  "BhS1", "BhS2",  "BhS3",  "BhS4",  "BhS5",  "BhS6",  "BhS7",  "BhS8", "BhS9","BhS10"))

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/apr21_subgspecTEs.tiff", units="in", width=6, height=3.5, res=300)
ggplot(bychr_totals, aes(x=chrom, y=perc, fill=genome))+
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=16),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    legend.position = "top",
    axis.title.x = element_text(size=17),
    axis.title.y = element_text(size=17, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=14, angle = -90, hjust=0.05, vjust=0.5),
    axis.text.y = element_text(size=14),
    panel.background = element_blank(),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
    axis.line = element_line(colour = "gray"),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))+
    ylab("Percent of TEs that are\nsubgenome-specific") +
    xlab("Chromosome") +
    scale_fill_manual(
      name = "Genome",
      values = c("#8ED6EE", "#574C47"),
      labels=c("ABR113", "Bhyb26"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(bychr_totals$perc)*1.1), labels = scales::percent_format(accuracy = 5L))
dev.off()

t.test(
  subset(bychr_totals, genome=="ABR113")$perc, 
  subset(bychr_totals, genome=="Bhyb26")$perc, 
  paired = TRUE, alternative = "two.sided")






