library(ggplot2)
library(tidyr)
library(stringr)


#manually copying output from analyze_counts_by_subg.py,
#see my notebook Nov.29 2021
countdat <- data.frame(
  'tissue' = c('Root', 'Leaf', 'Floret', 'Callus'),
  'D' = c(502979, 622342, 709650, 221374),
  'S' = c(485971, 571728, 819356, 205699)
)
countdat <- gather(countdat, subg, num_counts, D:S, factor_key=TRUE)





#NEW: red line indicates what percent of genes are from the BhS subgenome
#Bhyb26:
#wc -l BhDprimTrlengths.txt
#27265 BhDprimTrlengths.txt
# wc -l BhSprimTrlengths.txt
#26537 BhSprimTrlengths.txt
# 26537/(26537+27265) = 0.49323
tiff("/home/virginia/Documents/school/vogelLab/notebook/2021/nov30_subgcount.tiff", units="in", width=6, height=4, res=300)
ggplot(countdat, aes(x=tissue, y=num_counts, fill=subg))+
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_hline(yintercept = 0.49323, color='red', size=1)+ 
  theme(
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20, margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x = element_text(size=17),
    axis.text.y = element_text(size=17),
    panel.background = element_blank(),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
    axis.line = element_line(colour = "gray"),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))+
  ylab("Percent of total counts")+
  scale_fill_manual(
    name = "Read origin",
    values =c("#E69A8DFF", "#5F4B8BFF"), 
    labels=c("BhD transcript", "BhS transcript"))+
  scale_y_continuous(labels = scales::percent)
dev.off()


#manually copying output from horseraceUpdatedNEW.py,
#see my notebook Feb. 2 2022
horserace <- data.frame(
  'tissue' = c('Root', 'Leaf', 'Floret', 'Callus'),
  'favorsD' = c(6488, 5898, 6704, 5691),
  'total_pairs' = c(12937, 12010, 13376, 11414)
)

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb2_horserace.tiff", units="in", width=6.5, height=3, res=300)
#tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb9_ABR113_horserace.tiff", units="in", width=5, height=3, res=300)
#Bhyb26:
ggplot(horserace, aes(x=tissue, y=num_genes, fill=favored_subg))+
#ABR113:
#ggplot(horserace, aes(x=tissue, y=mean, fill=favored_subg))+
  geom_bar(stat = 'identity', position = 'dodge') +
  theme(
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=19, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=17),
    axis.text.y = element_text(size=15),
    panel.background = element_blank(),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
    axis.line = element_line(colour = "gray"),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))+
  ###FOR ABR113:
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))+
  ylab("Number of gene pairs") +
  scale_fill_manual(
    name = "Favored homeolog",
    values = c("sienna1", "skyblue1"),
    labels=c("BhD", "BhS"))
dev.off()



  
