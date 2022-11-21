suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

setwd('/home/virginia/Documents/school/vogelLab/notebook/2022')

TPMs <- read.table('Bhyb26_TPMs_w_pseudogenes.txt', sep='\t', header = TRUE, stringsAsFactors = FALSE)
conserved_df <- read.table('/home/virginia/Documents/school/vogelLab/notebook/2021/conservedgenes.txt', sep='\t', header = TRUE, stringsAsFactors = FALSE)
conserved_v <- c(conserved_df$Bhyb26Dgene, conserved_df$Bhyb26Sgene)
TPMs$Genetype <- "Other"
TPMs[which(TPMs$Gene %in% conserved_v),]$Genetype <-"Conserved genes"
TPMs[which(startsWith(TPMs$Gene, "Bradi") | startsWith(TPMs$Gene, "Brast")),]$Genetype <-"Putative pseudogenes"

mean(subset(TPMs, Genetype=="Conserved genes")$Total)
mean(subset(TPMs, Genetype=="Putative pseudogenes")$Total)

median(subset(TPMs, Genetype=="Conserved genes")$Total)
median(subset(TPMs, Genetype=="Putative pseudogenes")$Total)

#Not surprisingly, t-test is not significant.
t.test(subset(TPMs, Genetype=="Putative pseudogenes")$Total, subset(TPMs, Genetype=="Conserved genes")$Total)
#p-value = 0.7483

##### Later add-on: get mean/median TPM of ABR113 pseudogenes #####
ABR113_TPMs <- read.table('ABR113_TPMs_w_pseudogenes.txt', sep='\t', header = TRUE, stringsAsFactors = FALSE)
ABR113_pseudos <- subset(ABR113_TPMs, startsWith(ABR113_TPMs$Gene, "Bradi") | startsWith(ABR113_TPMs$Gene, "Brast"))
median(ABR113_pseudos$TPM)
#result: 0

#Get median TPM of Bhyb26 leaf data to make more comparable with ABR113
median(subset(TPMs, Genetype=="Conserved genes")$Leaf)
#result: 0.866
median(subset(TPMs, Genetype=="Putative pseudogenes")$Leaf)
#result: 0

mean(ABR113_pseudos$TPM)
# [1] 2.229461
mean(subset(TPMs, Genetype=="Conserved genes")$Leaf)
# [1] 18.09707
mean(subset(TPMs, Genetype=="Putative pseudogenes")$Leaf)
# [1] 10.61333

#This is why I originally didn't include ABR113--you can't compare TPMs across experiments...
#We can compare the ABR113 pseudogenes with the ABR113 conserved genes
gff_w_ogs <- read.table('/home/virginia/Documents/school/vogelLab/notebook/2021/gffWithOgs.txt', sep='\t', header = TRUE, stringsAsFactors = FALSE)
conserved_ABR113_gff <- ogs[which(gff_w_ogs$synOg %in% conserved_df$synOg & startsWith(gff_w_ogs$genome, "ABR113")),]
ABR113_TPMs$Genetype <- "Other"
ABR113_TPMs$Gene <- gsub(".1$", "", ABR113_TPMs$Gene) #remove the ".1" from end of geneIDs to be able to match to genespace gff
ABR113_TPMs[which(ABR113_TPMs$Gene %in% conserved_ABR113_gff$id),]$Genetype <-"Conserved genes"
ABR113_TPMs[which(startsWith(ABR113_TPMs$Gene, "Bradi") | startsWith(ABR113_TPMs$Gene, "Brast")),]$Genetype <-"Putative pseudogenes"
mean(subset(ABR113_TPMs, Genetype=="Conserved genes")$TPM)
# [1] 6.647757


#If we remove the single most highly expressed gene, does that make the means significantly different?
pseudoTPMs <- subset(TPMs, Genetype=="Putative pseudogenes")
pseudoTPMs <- pseudoTPMs[order(-pseudoTPMs$Total),] 
censored_pseudos <- pseudoTPMs[-1,]
# > mean(pseudoTPMs$Total)
# [1] 75.1521
# mean(censored_pseudos$Total)
# [1] 33.16571
res <- sapply(seq(length(j4)), function(n) mean(j4[[n]][,1]) < mean(censored_pseudos$Total) ) #this is the mean total TPM for the pseudogenes
sum(res)
t.test(censored_pseudos$Total, subset(TPMs, Genetype=="Conserved genes")$Total)
#p-value = 0.02205
#Apparently it does

p2 <- ggplot(subset(TPMs, Genetype!="Other"), aes(x=Total, fill=Genetype))+
  geom_density(color = "black", alpha = 0.5) + 
  scale_fill_manual( values = c("#E69F00", "#56B4E9"))+ #labels = c("Putative pseudogenes", "Control genes"),
  xlim(c(0, 23))+
  ylab("Density")+
  xlab("Transcripts Per Million (TPM)")+ #Note this is the total TPM from all 4 libraries
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=16),
    legend.position=c(0.7, 0.85),
    legend.key=element_blank(),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.box.background = element_rect(color="gray", size=1),
    legend.margin=margin(c(1,3,3,3)),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb18_TPMs.pseudos.tiff", units="in", width=5.5, height=3, res=300)
p2+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(ggplot_build(p2)$data[[1]]$ymax)*1.1))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 23*1.1))
dev.off()


###### Later add-on: include ABR113 ######

all_TPMs <- TPMs[,c("Gene", "Leaf", "Genetype")]
all_TPMs$Genome <- "Bhyb26"
names(all_TPMs) <- c("Gene", "TPM", "Genetype", "Genome")
ABR113_TPMs$Genome <- "ABR113"
all_TPMs <- rbind(all_TPMs, ABR113_TPMs)
all_TPMs$gg <- paste0(all_TPMs$Genome, " ", all_TPMs$Genetype)
all_TPMs <- subset(all_TPMs, Genetype!="Other")

p2 <- ggplot(subset(all_TPMs, Genome=="ABR113"), aes(x=TPM, fill=gg))+
  geom_density(color = "black", alpha = 0.5) + 
  scale_fill_manual( values = c("#E69F00", "#56B4E9"))+ #labels = c("Putative pseudogenes", "Control genes"),
  xlim(c(0, 4))+
  ylab("Density")+
  xlab("Transcripts Per Million (TPM)")+ 
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=16),
    legend.position=c(0.63, 0.85),
    legend.key=element_blank(),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.box.background = element_rect(color="gray", size=1),
    legend.margin=margin(c(1,3,3,3)),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/apr18_TPMs.pseudos.ABR113.tiff", units="in", width=5.5, height=3, res=300)
p2+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(ggplot_build(p2)$data[[1]]$ymax)*1.1))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4*1.1))
dev.off()


p3 <- ggplot(subset(all_TPMs, Genome=="Bhyb26"), aes(x=TPM, fill=gg))+
  geom_density(color = "black", alpha = 0.5) + 
  scale_fill_manual( values = c("#E69F00", "#56B4E9"))+ #labels = c("Putative pseudogenes", "Control genes"),
  xlim(c(0, 4))+
  ylab("Density")+
  xlab("Transcripts Per Million (TPM)")+ 
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=16),
    legend.position=c(0.63, 0.85),
    legend.key=element_blank(),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.box.background = element_rect(color="gray", size=1),
    legend.margin=margin(c(1,3,3,3)),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/apr18_TPMs.pseudos.Bhyb26.tiff", units="in", width=5.5, height=3, res=300)
p3+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(ggplot_build(p3)$data[[1]]$ymax)*1.1))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4*1.1))
dev.off()



















