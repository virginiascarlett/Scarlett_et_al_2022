library(ggplot2)
library(tidyr)
options(scipen=999) #to prevent scientific notation

get_initial_data <- function(genome) {
  myages <- read.csv(paste0('finalAges.txt_', genome), sep='\t', header=TRUE, stringsAsFactors = FALSE) #produced in step 7 of my TE pipeline (EMBOSS)
  myfrags <- read.csv(paste0('allfragments.classified_', genome), sep='\t', header=TRUE, stringsAsFactors = FALSE)
  myages$Classification <- factor(myfrags[match(myages$TEcopy, myfrags$TEcopycode),]$classification)
  myages$Genome <- genome
  myages$Exemplar <- factor(myfrags[match(myages$TEcopy, myfrags$TEcopycode),]$exemplar)
  return(myages)
}
get_lifespan_df <- function(mydata, mygenome) {
  famsizes <- table(mydata$Exemplar)
  t <- famsizes[which(famsizes >= 2)] #we can only get lifespans for families that have at least 2 members, 
  #which is about one-third of the total number of families!
  mins <- sapply(names(t), function(name) min(subset(mydata, Exemplar==name)$InsertionTime))
  maxs <- sapply(names(t), function(name) max(subset(mydata, Exemplar==name)$InsertionTime))
  df <- data.frame(Exemplar=names(t), Youngest=mins, Oldest=maxs, Genome=mygenome)
  df[,c(2:3)] <- df[,c(2:3)]/1000000
  df <- df[order(df$Oldest, decreasing = T),]
  return(df)
}

setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')
ABR113 <- get_initial_data('ABR113')
Bhyb26 <- get_initial_data('Bhyb26')
ABR114 <- get_initial_data('ABR114')
Bd21 <- get_initial_data('Bd21')
Bd1_1 <- get_initial_data('Bd1_1')

#% of dateable LTR-RTs that are less than 1.4 MY old
nrow(subset(ABR113, InsertionTime < 1400000))/nrow(ABR113)
nrow(subset(Bhyb26, InsertionTime < 1400000))/nrow(Bhyb26)


final <- rbind(ABR113, Bhyb26)
final$famgen <-  paste0(final$Classification, '_', final$Genome)

ABR113_LS <- get_lifespan_df(ABR113, "ABR113")
Bhyb26_LS <- get_lifespan_df(Bhyb26, "Bhyb26")
Bd1_1_LS <- get_lifespan_df(Bd1_1, "Bd1_1")
Bd21_LS <- get_lifespan_df(Bd21, "Bd21")
ABR114_LS <- get_lifespan_df(ABR114, "ABR114")

# Our lifespans dataframe has one row per family (per genome), with the min and max age of elements in that family.
# These will appear as vertical lines on the plot.
lifespans_all <- rbind( ABR113_LS, Bhyb26_LS, Bd1_1_LS, Bd21_LS, ABR114_LS )
lifespans_all$Intercept <- seq(0, 20, length.out=nrow(lifespans_all))
# Our allages_all dataframe has one row per element (per genome), so that ALL TEs will appear as a point on the plot
allages_all <- rbind(merge(lifespans_all, Bd21),
                     merge(lifespans_all, Bd1_1),
                     merge(lifespans_all, ABR114),
                     merge(lifespans_all, Bhyb26),
                     merge(lifespans_all, ABR113))

#allages_all <- rbind(merge(lifespans_all, Bhyb26), merge(lifespans_all, ABR113))


tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/mar28_lifespans.tiff", units="in", width=7, height=4.25, res=300)
ggplot()+
  geom_segment(data=lifespans_all, aes(x = Intercept, y = Youngest, xend = Intercept, yend = Oldest, color=Genome))+
  scale_color_manual(labels = c("ABR113", "Bhyb26", "Bd1-1", "Bd21", "ABR114"), values = c("deepskyblue4", "seagreen3", "plum2", "khaki3", "purple3"))+
  geom_hline(aes(yintercept=1.4), size=0.65, color="seagreen3")+
  geom_hline(aes(yintercept=0.14), size=0.65, color="deepskyblue4")+ 
  geom_point(data=allages_all, aes(x=Intercept, y=InsertionTime/1000000, color=Genome), size=1)+
  ylab("Insertion Time (MYA)")+
  xlab("Intact LTR-RTs by Family")+
  theme(
    legend.title=element_blank(), 
    legend.text=element_text(size=14),
    legend.position = "top",
    #legend.position=c(0.95, 0.9),
    legend.key=element_blank(),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.box.background = element_rect(color="gray"),
    legend.margin=margin(c(1,3,3,3)),
    axis.title.x = element_text(size=17),
    axis.title.y = element_text(size=17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
    axis.line = element_line(colour = "gray"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  scale_y_continuous(expand = c(0, 0))+ #limits = c(0, 10),
  coord_cartesian(ylim=c(0, 10), xlim=c(0, 20))+
  guides(color = guide_legend(override.aes = list(size = 1)))
dev.off()
