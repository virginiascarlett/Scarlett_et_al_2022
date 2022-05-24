library(ggplot2)

setwd("/home/virginia/Documents/school/vogelLab/notebook/2021/lost_genes")

#Bhyb26 conserved genes versus diploid orthologs
mydata <- read.table('dNdSbyGene.all_conserved_genes.txt', sep='\t', header = FALSE, stringsAsFactors = FALSE)
Bhyb26allcons <- mydata[!mydata$V2==99,] #99 indicates there were no synonymous substitutions

#Bhyb26 lost genes versus diploid orthologs
mydata <- read.table('dNdSbyGene.dipsvsBhyb26.txt', sep='\t', header = FALSE, stringsAsFactors = FALSE)
Bhyb26vsdip_lost <- mydata[!mydata$V2==99,] #99 indicates there were no synonymous substitutions

# Now let's do ABR113
mydata <- read.table('/home/virginia/Documents/school/vogelLab/notebook/2021/lost_genes/dNdSbyGene.dipsvsABR113.txt', sep='\t', header = FALSE, stringsAsFactors = FALSE)
ABR113vsdip_lost <- mydata[!mydata$V2==99,]

final <- data.frame(value=Bhyb26allcons$V2, group="Conserved genes")
final <- rbind(final, data.frame(value=ABR113vsdip_lost$V2, group="ABR113 lost genes"))
final <- rbind(final, data.frame(value=Bhyb26vsdip_lost$V2, group="Bhyb26 lost genes"))

                              
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/mar30_dNdS.tiff", units="in", width=4, height=3, res=300)
ggplot(final, aes(x=value, fill=group, linetype=group) )+
  #geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+
  geom_density(color = "black", alpha = 0.3) + 
  scale_fill_manual(name="my legend", labels = c("Conserved genes, Bhyb26", "Putative pseudogenes, ABR113", "Putative pseudogenes, Bhyb26"), values = c( "NA", "#8AB8C5", "#3E282B"))+ #c( "#E69F00", "NA", "#56B4E9")
  scale_linetype_manual(name="my legend", labels = c("Conserved genes, Bhyb26", "Putative pseudogenes, ABR113", "Putative pseudogenes, Bhyb26"), values=c("dashed", "solid", "solid"))+
  xlim(c(0, 2))+
  ylab("Density")+
  xlab("dN/dS")+
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=10.5),
    legend.position=c(0.635, 0.865),
    legend.key=element_blank(),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.box.background = element_rect(color="gray", size=1),
    legend.margin=margin(t=-1, unit='mm'), #margin(t=-2.1, unit='mm'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

dev.off()


