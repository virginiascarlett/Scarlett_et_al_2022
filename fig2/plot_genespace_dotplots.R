suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
})
options(scipen=999) #to prevent scientific notation

#This script is broken up into different sections
#Run the section that corresponds to the type of plot you want to make
#(hybridum vs. hybridum, BhS vs. stacei, etc.)

#We are just plotting the anchor genes, or synteny constrained orthologs

roundUp <- function(x,m) m*ceiling(x / m)

get_tp <- function(possiblefname1, possiblefname2) {
  return(
    if (file.exists(possiblefname1)) {
    d <- fread(possiblefname1, na.strings = c("", "NA"))
    tp <- subset(d, isAnchor)
  } else {
    d <- fread(possiblefname2, na.strings = c("", "NA"))
    tp <- subset(d, isAnchor)
    tp <- data.frame(
      gen1=tp$gen2,
      gen2=tp$gen1,
      start1=tp$start2,
      start2=tp$start1,
      end1=tp$end2,
      end2=tp$end1,
      chr1=tp$chr2,
      chr2=tp$chr1,
      blkID=tp$blkID
    )
  }
  )
}

plot_genome <- function(plottingdf, xvar="start1", yvar="start2", xtitle, ytitle) {
  return(
    ggplot(plottingdf, aes_string(x = xvar, y = yvar))+
    #geom_point(aes(color=blkID), pch = ".", size=2)+
      geom_point(pch = ".", size=2, color="blue")+
    theme(
      legend.position = "none",
      axis.title.x = ggtext::element_markdown(size=14),
      axis.title.y = ggtext::element_markdown(size=14),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.border=element_rect(color="gray", fill=NA),
      #panel.background = element_blank(),
      panel.background = element_rect(color="gray"),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(color="white", fill="white"),
      strip.text.y.left = element_text(size = 12, angle=0),
      strip.text.x = element_text(size = 12, angle=-70),
      axis.ticks = element_blank(),
      plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
    scale_x_continuous(name=xtitle)+
    scale_y_continuous(name=ytitle)+
    facet_grid(chr2 ~ chr1, scale = "free", space = "free", as.table = F, switch="both")
  )
}

plot_chromosome <- function(plottingdf, xvar="start1", yvar="start2", xtitle, ytitle, plottingchr1max, plottingchr2max) {
  return(
    ggplot(plottingdf, aes(x = start1, y = start2))+
      geom_point(aes(color=blkID), pch = ".")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.border=element_rect(color="gray", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
        panel.spacing = unit(0, "lines"),
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
      facet_grid(chr2 ~ chr1, scale = "free", space = "free", as.table = F, switch="both")+
      scale_x_continuous(name=xtitle, breaks = seq(0, plottingchr1max, by=2000000), labels=seq(0, plottingchr1max, by=2000000)/1000000, expand = c(0, 0))+
      scale_y_continuous(name=ytitle, breaks = seq(0, plottingchr2max, by=2000000), labels=seq(0, plottingchr2max, by=2000000)/1000000, expand = c(0, 0))+
      coord_cartesian(xlim = c(0, plottingchr1max), ylim = c(0, plottingchr2max) )
  )
}

plot_chromosomeNEW <- function(plottingdf, xvar="start1", yvar="start2", xtitle, ytitle, plottingchr1max, plottingchr2max) {
  return(
    ggplot(plottingdf, aes(x = start1, y = start2))+
      geom_point(size=0.5, color="blue")+
      theme(
        legend.position = "none",
        axis.title.x = ggtext::element_markdown(size=28),
        #axis.title.y = element_text(size=28),
        axis.title.y = ggtext::element_markdown(size=28),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border=element_rect(color="gray10", fill=NA),
        panel.background = element_rect(fill="gray90"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"),
        #panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
        panel.spacing = unit(0, "lines"),
        #strip.background=element_blank(),
        #strip.text=element_blank(),
        strip.background = element_rect(color="white", fill="white"),
        strip.text.y.left = element_text(size = 24, angle=0),
        strip.text.x = element_text(size = 24),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
      facet_grid(chr2 ~ chr1, scale = "free", space = "free", as.table = F, switch="both")+
      scale_x_continuous(name=xtitle, breaks = seq(0, plottingchr1max, by=2000000), labels=seq(0, plottingchr1max, by=2000000)/1000000, expand = c(0, 0))+
      scale_y_continuous(name=ytitle, breaks = seq(0, plottingchr2max, by=2000000), labels=seq(0, plottingchr2max, by=2000000)/1000000, expand = c(0, 0))+
      coord_cartesian(xlim = c(0, plottingchr1max), ylim = c(0, plottingchr2max) )
  )
}






################### MAKE DOT PLOT S SUBGENOME OR A SPECIFIC BhS CHROMOSOME VERSUS B. STACEI ################### 

hybgenome = 'Bhyb26S' #EDIT THIS LINE

tp <- get_tp(
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/Bstacei_", hybgenome, "_synHits.txt.gz"),
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/", hybgenome, "_Bstacei_synHits.txt.gz")
)

tp$chr2 <- factor(tp$chr2, levels = c("BhS1","BhS2","BhS3","BhS4","BhS5","BhS6","BhS7","BhS8","BhS9","BhS10"))

#p <- plot_genome(plottingdf=tp, xtitle=expression(paste(italic("B. stacei"), " ABR114")), ytitle=paste0("*B. hybridum* ", hybgenome)) 
p <- plot_genome(plottingdf=tp, xtitle="*B. stacei* ABR114", ytitle=paste0("*B. hybridum* ", hybgenome)) 
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/Bstacei_", hybgenome, "_dotplot.tiff"), plot=p, units="in", width=5, height=5, dpi=300)


mychrom <- "BhS7"
toPlot <- subset(tp, chr2 == mychrom)
chr1max <- roundUp(max(toPlot$end1), 1000000)
chr2max <- roundUp(max(toPlot$end2), 1000000)
p <- plot_chromosome(plottingdf=toPlot, xtitle=paste0("ABR114 ", gsub('BhS', 'Bs', mychrom), " (Mb)"), ytitle=paste0(hybgenome, " ", mychrom, " (Mb)"), plottingchr1max=chr1max, plottingchr2max=chr2max)
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/dotplots_by_chromosome/", mychrom, "_Bstacei_", hybgenome, ".tiff"), plot=p, units="in", width=5, height=5, dpi=300)

#SPECIAL ADD-ON: TWO CHROMOSOMES 
toPlot <- subset(tp, chr2 == 'BhS7' | chr2 == 'BhS8')
chr1max <- roundUp(max(toPlot$end1), 1000000)
chr2max <- roundUp(max(toPlot$end2), 1000000)
#AXIS TITLES ARE A WORK IN PROGRESS
p <- plot_chromosomeNEW(plottingdf=toPlot, xtitle="*B. stacei* ABR114", ytitle=paste0("*B. hybridum* ", hybgenome), plottingchr1max=chr1max, plottingchr2max=chr2max)
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/dotplots_by_chromosome/7and8_Bstacei_", hybgenome, ".tiff"), plot=p, units="in", width=5, height=5, dpi=300)












################### NOW LET'S DO THE SAME FOR D GENOMES ################### 

hybgenome = 'Bhyb26D' #EDIT THIS LINE

tp <- get_tp(
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/Bdistachyon_", hybgenome, "_synHits.txt.gz"),
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/", hybgenome, "_Bdistachyon_synHits.txt.gz")
)

tp$chr2 <- factor(tp$chr2, levels = c("BhD1","BhD2","BhD3","BhD4","BhD5"))


p <- plot_genome(plottingdf=tp, xtitle="*B. distachyon* Bd21", ytitle=paste0("*B. hybridum* ", hybgenome)) 
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/Bdistachyon_", hybgenome, "_dotplot.tiff"), plot=p, units="in", width=5, height=5, dpi=300)


mychrom <- "BhD5"
toPlot <- subset(tp, chr2 == mychrom)
chr1max <- roundUp(max(toPlot$end1), 1000000)
chr2max <- roundUp(max(toPlot$end2), 1000000)
p <- plot_chromosome(plottingdf=toPlot, xtitle=paste0('Bd21 ', gsub('BhD', 'Bd', mychrom), " (Mb)"), ytitle=paste0(hybgenome, " ", mychrom, " (Mb)"), plottingchr1max=chr1max, plottingchr2max=chr2max)
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/dotplots_by_chromosome/", mychrom, "_Bdistachyon_", hybgenome, ".tiff"), plot=p, units="in", width=5, height=5, dpi=300)











  


################### WHY NOT HYBRIDUM VS HYBRIDUM TOO? JUST BhS FOR NOW ################### 

mygenome1 = 'ABR113S' #EDIT THIS LINE
mygenome2 = 'Bd28S' #EDIT THIS LINE

tp <- get_tp(
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/", mygenome1, "_", mygenome2, "_synHits.txt.gz"),
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/results/", mygenome2, "_", mygenome1, "_synHits.txt.gz")
)

#IF YOU ARE LOOKING AT S SUBGENOMES:
tp$chr1 <- factor(tp$chr2, levels = c("BhS1","BhS2","BhS3","BhS4","BhS5","BhS6","BhS7","BhS8","BhS9","BhS10"))
tp$chr2 <- factor(tp$chr2, levels = c("BhS1","BhS2","BhS3","BhS4","BhS5","BhS6","BhS7","BhS8","BhS9","BhS10"))


p <- plot_genome(plottingdf=tp, xtitle=paste0("B. hybridum ", mygenome1), ytitle=paste0("B. hybridum ", mygenome2)) 
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/", mygenome1, "_", mygenome2, "_dotplot.tiff"), plot=p, units="in", width=5, height=5, dpi=300)

  mychrom <- "BhS3"
  toPlot <- subset(tp, chr2 == mychrom)
  chr1max <- roundUp(max(toPlot$end1), 1000000)
  chr2max <- roundUp(max(toPlot$end2), 1000000)
  p <- plot_chromosome(plottingdf=toPlot, xtitle=paste0(mygenome1, " ", mychrom, " (Mb)"), ytitle=paste0(mygenome2, " ", mychrom, " (Mb)"), plottingchr1max=chr1max, plottingchr2max=chr2max)
  ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/dotplots_by_chromosome/", mychrom, "_", mygenome1, "_", mygenome2, ".tiff"), plot=p, units="in", width=5, height=5, dpi=300)














################### WE KNEW THIS WAS COMING: DISTACHYON VS DISTACHYON ################### 


mygenome2 = 'Bd301' #EDIT THIS LINE. Bd301 or Bd11

tp <- get_tp(
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_distachyon_new/results/Bdistachyon_", mygenome2, "_synHits.txt.gz"),
  paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_distachyon_new/results/", mygenome2, "_Bdistachyon_synHits.txt.gz")
)

#VERY IMPORTANT!! YOU MUST MANUALLY INPUT THE AXIS TITLE: ADD THE DASH IN Bd1-1 or Bd30-1!!
p <- plot_genome(plottingdf=tp, xtitle="B. distachyon Bd21", ytitle="B. distachyon Bd30-1")
ggsave(paste0("/home/virginia/Documents/school/vogelLab/GENESPACE_distachyon_new/Bdistachyon", "_", mygenome2, "_dotplot.tiff"), plot=p, units="in", width=5, height=5, dpi=300)


























