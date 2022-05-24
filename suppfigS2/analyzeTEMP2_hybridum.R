#A script to plot temp2 output produced by the mcclintock pipeline.
#We are trying to get a dataframe with one column for the sample id,
#another column for the TE name, and another column for how many good insertions that TE has in that sample
#(how many times it appears in the filtered df).

#Note: TEMP2 provides a comma-separated list  of TEs if the read matches multiple TEs. 
#I am choosing to keep the longest hit. I don't really see a better way of doing this
#without seriously complicating the workflow. It's simply difficult to summarize the
#results when most of the results are ambiguous.

#I also use this script to check whether library depth or quality correlates with total # of TIPs.
#To do this part, need to run FASTQC on all samples, grab the file fastqc_data.txt from each,
#and then run my script TEMP2QC.R before running this script.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyr)
  library(corrplot)
  library(stringr)
})
options(scipen=999) #to prevent scientific notation

nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
mycolors2 <- colorRampPalette(brewer.pal(9, "Set1"))(11)
#test <- colorRampPalette(c("tomato3", "gold", "olivedrab1", "steelblue1", "saddlebrown"))(17)
test <- colorRampPalette(c("tomato3", "deepskyblue4", "seagreen3", "plum2", "gold", "purple3"))(17)
test2 <- c(colorRampPalette(c("tomato3", "gold", "seagreen3", "deepskyblue4","purple3", "plum2"))(10), "gray50")

read_my_file <- function(filename) {
  d <- read.table(filename, sep='\t', header = FALSE, stringsAsFactors = FALSE)#  
  return(d)
}
filter_my_df <- function(df){
  return( subset(df, df$V7=='1p1' & df$V5 > 0.20) )
}
get_TIP_contribution <- function(genome, TEname) { #what fraction of total TIPs in this genome are from this TE family?
  thisfam <- subset(final_data, TE==TEname & Sample==genome)$Freq
  allfams <- sum(subset(final_data, Sample==genome)$Freq)
  #final_data[which(final_data$TE==TEname),]$TIPcontrib <<- thisfam/allfams
  return(thisfam/allfams)
}
get_longest <- function(s) {
  myTElist <- unlist(strsplit(as.character(s), ","))
  myTElist <- unlist(unname(sapply(myTElist, function(myTE) unlist(strsplit(myTE, ":"))[1] )))
  lens <- len_dat[match(myTElist, len_dat$V1),]
  longest <- lens[which(lens$V2==max(lens$V2)), "V1"]
  if (length(longest > 1)) {return(longest[1])
  } else {return(longest)}
}

setwd("/home/virginia/Documents/school/vogelLab/notebook/2021/TEpolymorphisms/individual_annots") #Set to directory with files
files <- list.files(pattern = "\\.bed") #Collects all .bed files in directory

list_of_dataframes <- lapply(files, function(n) read_my_file(n))
#If you are getting the error "no lines available in input" then you have one or more empty files.
#You can try mytest <- lapply(files, function(n) tryCatch(read_my_file(n), error=function(e) NULL
#And then do sapply(mytest, function(x) x[1,c(1:3)]) to see which files are NULL
#for each dataframe, create a column for the sample ID
sampleIDs <- unlist(unname(sapply(files, function(filename) gsub("_R1.insertion.bed", "", gsub("Bhybridum", "", filename)))))
junk <- lapply(seq(length(files)), function(n) list_of_dataframes[[n]]$Sample <<- sampleIDs[n])
my_data <- rbindlist(list_of_dataframes) 
my_data$V4 <- unlist(lapply(my_data$V4, function(data) gsub("UNDERSCORE", "_", gsub("DASH", "-", data))))
filtered_data <- filter_my_df(my_data)

#PULL THE LONGEST EXEMPLAR 
len_dat <- read.csv('exemplarlengths.ABR113.txt', sep='\t', header = FALSE, stringsAsFactors = FALSE)
len_dat$V1 <- sapply(len_dat$V1, function(s) unlist(strsplit(s, "#"))[1] )
picks <- unlist(unname(sapply(filtered_data$V4, function(x) get_longest(x))))
filtered_data$TE <- picks
filtered_data <- filtered_data[,c("Sample","TE")]

final_data <- as.data.frame(table(filtered_data))
#For master annotation:
# samplekey <- read.csv('README.concise2', header = FALSE, stringsAsFactors = FALSE)
# samplekey$V1 <- unlist(unname(sapply(samplekey$V1, function(x) strsplit(x, "\\.")[[1]][1])))
#For invidiual annotations:
samplekey <- read.csv('TEMP2samplekey.txt', header = FALSE, stringsAsFactors = FALSE)

final_data$Sample <-samplekey[match(final_data$Sample, samplekey$V1),]$V2
final_data$TIPcontrib <- unname(mapply(function(g, t) get_TIP_contribution(g, t), final_data$Sample, final_data$TE))
accns_to_keep <- c(samplekey$V2[-which(samplekey$V2=="Bd3-1xSta5_S4" | samplekey$V2=="Bd3-1_Sta5_concat")])
final_data <- subset(final_data, Sample %in% accns_to_keep)


##### QC check #######
get_total_TIPs <- function(mysample){
  mini <- subset(final_data, Sample==mysample)
  return(sum(mini$Freq))
}
totals <- sapply(unique(final_data$Sample), function(s) get_total_TIPs(s))
QC <- read.csv(file="/home/virginia/Documents/school/vogelLab/notebook/2021/TEpolymorphisms/QC/QCinfo.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
get_ID <- function(mystr) {
  s <- substr(mystr, 19, nchar(mystr))
  return(substr(s, 1, nchar(s)-4))
}
rownames(QC) <- paste0("mcclintock", unlist(unname(sapply(rownames(QC), function(mystr) get_ID(mystr)))))
rownames(QC) <- samplekey[match(rownames(QC), samplekey$V1),]$V2
QC$total_TIPs <- totals[match(rownames(QC), names(totals))]
QC <- QC[!is.na(QC$total_TIPs),]
ggplot(QC, aes(x = total_TIPs, y=medianQ)) + 
  geom_point()+
  geom_smooth(method=lm)
ggplot(QC, aes(x = total_TIPs, y=depth)) + 
  geom_point()+
  geom_smooth(method=lm)

m <- lm(medianQ ~ total_TIPs, QC)
summary(m)
########################

#Master:
#The top 39% of total TIPs come from 10 TE families
#The top 75% of total TIPs come from 100 TE families
#Individual:
#The top 52.5% of total TIPs come from 10 TE families
#The top 87.0%% of total TIPs come from 100 TE families
famtotals <- aggregate(Freq~TE, final_data, sum)
famtotals <- famtotals[order(famtotals$Freq, decreasing = TRUE),]
sum(famtotals[1:10,'Freq'])/sum(famtotals$Freq)
sum(famtotals[1:100,'Freq'])/sum(famtotals$Freq)

#Master:
#No single TE family contributes more than 16% of a genome's total TIPs
#Individual:
#No single TE family contributes more than 25.4% of a genome's total TIPs
fams_by_contrib <- final_data[order(final_data$TIPcontrib, decreasing = TRUE),]
  
classkey <- read.csv('/home/virginia/Documents/school/vogelLab/scripts/TEs/three_letter_codes.txt', header=FALSE, skip=1, stringsAsFactors = FALSE)
TElib <- read.csv('/home/virginia/Documents/school/vogelLab/notebook/2021/TElib.ABR113.fa', sep='\n', header=FALSE, stringsAsFactors = F)
TElib <- TElib$V1[c(TRUE, FALSE)] #get every other element, namely the library headers
temp <- unlist(strsplit(TElib, "#"))
TEdf <- data.frame(
  Family = gsub(">", "", temp[c(TRUE, FALSE)]),
  Libclass = temp[c(FALSE, TRUE)],
  stringsAsFactors = FALSE
)
TEdf$Libclass <- sapply(strsplit(TEdf$Libclass, " "), `[[`, 1) #get everything before the space, e.g. for DNA/MULE-MuDR RepbaseID: MuDR-18_SBiXX just get DNA/MULE-MuDR
TEdf$Classification <- classkey[match(TEdf$Libclass, classkey$V1),]$V2
final_data$Classification <- TEdf[match(final_data$TE, TEdf$Family),"Classification"]

#Get the abundance of TE classes in the reference genome
fragdat <- read.csv('/home/virginia/Documents/school/vogelLab/notebook/2021/allfragments.classified_ABR113', sep='\t', header=TRUE, stringsAsFactors = FALSE)
genomeTEs <- fragdat[!duplicated(fragdat[,c('TEcopycode')]),c('TEcopycode', 'exemplar', 'classification')]
genomeTEs$fam <- unlist(unname(sapply(genomeTEs$exemplar, function(s) unlist(unname(strsplit(s, '#')))[1] )))
genomeTEs$superfam <- unlist(unname(sapply(genomeTEs$exemplar, function(s) unlist(unname(strsplit(unlist(unname(strsplit(s, '#')))[2], "\\s+")))[1] )))
Reference <- table(genomeTEs$classification)[order(unlist(table(genomeTEs$classification)))]
tally <- aggregate(Freq~Classification, final_data, sum)
TIPs <- tally$Freq
names(TIPs) <- tally$Classification
res <- as.data.frame(cbind(Reference, TIPs = TIPs[names(Reference)]))
res$Classification <- rownames(res)
res[is.na(res)] <- 0

#Plot the relative contribution of each TE class to all ABR113 TEs vs. 
#the relative contribution of each TE class to all TIPs
#(Also plotted with the relative contribution of each TE class to the 100 families with the most total TIPs 
#and the plot was almost exactly the same)
row.names(res) <- NULL
toPlot <- gather(res, Source, Total, c('Reference', 'TIPs'), factor_key=TRUE)
toPlot$Classification <- gsub("DNA_transposon", "DNA", toPlot$Classification)
toPlot$Classification <- gsub("LTR_retrotransposon_", "LTR_retro_", toPlot$Classification)
toPlot$Classification <- gsub("non_LTR_retrotransposon", "non_LTR_retro", toPlot$Classification)
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/mar28_hybridumTIPs1b.tiff", units="in", width=6.5, height=4.25, res=300)
ggplot(toPlot, aes(x = Source, y = Total, fill = Classification)) + 
  geom_bar(stat = "identity", position="fill")+ 
           #width= 0.6, #how wide (0-1) the bars are
           #color = "black", # the outline of the bars
          # size= 0.5, # the thickness of the outline
           #alpha=0.7) +
  #ylab("TIPs")+
  theme(
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size=13),
    legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=16), #angle=-70),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = test)+
  ylab("Percent of Total")+
  scale_x_discrete(labels=c("Reference" = "Reference TEs", "TIPs"="TIPs")) #"TIPs"="Transposon Insertion\nPolymorphisms (TIPs)"))

dev.off()

# tiff("/home/virginia/Documents/school/vogelLab/notebook/2021/nov3_hybridumTIPs.tiff", units="in", width=6.75, height=4.25, res=300)
# ggplot(toPlot, aes(x = Source, y = Total, fill = Classification)) + 
#   geom_bar(stat = "identity", position="fill")+ 
#   #width= 0.6, #how wide (0-1) the bars are
#   #color = "black", # the outline of the bars
#   # size= 0.5, # the thickness of the outline
#   #alpha=0.7) +
#   #ylab("TIPs")+
#   theme(
#     legend.position = "right",
#     legend.title = element_blank(), 
#     legend.text = element_text(size=13),
#     legend.spacing.y = unit(0, 'cm'),
#     legend.key.height = unit(0.2, 'cm'),
#     legend.key.width = unit(0.5, 'cm'),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size=16),
#     axis.text.x = element_text(size=16), #angle=-70),
#     axis.text.y = element_text(size=15),
#     panel.border=element_rect(colour="gray", fill=NA),
#     panel.background = element_blank(),
#     panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
#     panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
#     plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
#   scale_y_continuous(labels = scales::percent)+
#   scale_fill_manual(values = test)+
#   ylab("Percent of Total")+
#   scale_x_discrete(labels=c("Reference" = "Reference\nTEs", "TIPs"="Transposon Insertion\nPolymorphisms (TIPs)"))
# dev.off()

#Next, let's work on plotting the TE families with the highest median
#Next, create a vector of TE names in the order we want them (for now we will order by median frequency)
get_fam_median <- function(famname, df) {
  mini <- subset(df, TE==famname)
  return(median(mini$Freq))
}

famstats <- sapply(unique(final_data$TE), function(name) get_fam_median(name, final_data))
names(famstats) <- unique(final_data$TE)
famsordered <- names(famstats[order(unlist(famstats), decreasing=TRUE)])
topTEs <- head(famsordered, 15) #Display the top n TE families (those with the highest median TIPs)
toPlot <- final_data
toPlot$TE <- as.character(toPlot$TE)
toPlot[which(!(toPlot$TE %in% topTEs)),"TE"] <- "Other"
toPlot <- aggregate(toPlot['Freq'], by=toPlot[,c("Sample", "TE")], sum)

#Sloppy coding, just putting our favorite accessions first, one at a time
#MASTER:
#accnsordered <- c("Bhyb118-5", unique(toPlot$Sample)[-which(unique(toPlot$Sample)=="Bhyb118-5")])
#INDIVIDUAL:
accnsordered <- c("Bhyb118_5", unique(toPlot$Sample)[-which(unique(toPlot$Sample)=="Bhyb118_5")])
accnsordered <- c("Bhyb26", accnsordered[-which(accnsordered=="Bhyb26")])
accnsordered <- c("ABR113", accnsordered[-which(accnsordered=="ABR113")])
#omit mislabeled distachyon samples
accnsordered <- accnsordered[-which(accnsordered=="Bhyb123_1")]
accnsordered <- accnsordered[-which(accnsordered=="Bhyb130_1")]
accnsordered <- accnsordered[-which(accnsordered=="Bhyb162_1")]

#Plot by TIP contribution: Choose the top 10 TEs that are 
#the biggest contributors to that accession's TIPs in any line.
topTIPcontrs <- as.character(unique(fams_by_contrib$TE)[1:10])
toPlotContr <- final_data
toPlotContr$TE <- as.character(toPlotContr$TE)
toPlotContr[which(!(toPlotContr$TE %in% topTIPcontrs)),"TE"] <- "Other"
toPlotContr <- aggregate(toPlotContr['Freq'], by=toPlotContr[,c("Sample", "TE")], sum)
toPlotContr <- subset(toPlotContr, Sample %in% accnsordered)
toPlotContr$Sample <- factor(toPlotContr$Sample, levels=accnsordered)
toPlotContr$TE <- factor(toPlotContr$TE, levels=c(topTIPcontrs, "Other"))

fullname <- paste0(TEdf[match(toPlotContr$TE, TEdf$Family),]$Family, "#", TEdf[match(toPlotContr$TE, TEdf$Family),]$Libclass)
fullname <- gsub("Gypsy", "RLG", fullname)
fullname <- gsub("Copia", "RLC", fullname)
fullname <- gsub("NA#NA", "All other families", fullname)
toPlotContr$fullname <- factor(fullname, levels=
                                 c(unique(fullname)[-which(unique(fullname)=="All other families")], "All other families")
                               )

#tiff("/home/virginia/Documents/school/vogelLab/notebook/2021/nov2_hybridumTIPs2.tiff", units="in", width=8, height=4.25, res=300)
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/mar28_hybridumTIPs.tiff", units="in", width=8.5, height=4.25, res=300)
ggplot(toPlotContr, aes(x = Sample, y = Freq, fill = fullname)) + 
    geom_bar(stat = "identity", 
             width= 0.6, #how wide (0-1) the bars are
             color = "grey20", # the outline of the bars
             size= 0.5, # the thickness of the outline
             alpha=0.7) +
    ylab("Transposon Insertion Polymorphisms\n(TIPs)")+
    theme(
      #legend.position = "none",
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size=11),
      legend.spacing.y = unit(0.5, 'cm'),
      legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(0.5, 'cm'),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=14),
      axis.text.x = element_text(angle = -90, size=14, hjust=0.05, vjust=0.5),
      axis.text.y = element_text(size=14),
      panel.border=element_rect(colour="gray", fill=NA),
      panel.background = element_blank(),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
      panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
      plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    scale_fill_manual(values = test2) +
    scale_y_continuous(limits = c(0,600), expand = c(0, 0))
  dev.off()

# SCRAPS
#chisq <- chisq.test(res[,c(1,2)], simulate.p.value = T)
#round(chisq$residuals, 3)
#p-value = 0.0004998
#we also see that RLC TIPs have a high positive residual
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
#https://scholarworks.umass.edu/cgi/viewcontent.cgi?article=1269&context=pare


#fams_filt <- subset(famtotals, Freq >= 5)$TE
#final_small <- subset(final_data, TE %in% fams_filt)
#final_small$Classification <- genomeTEs[match(final_small$TE, genomeTEs$fam),"classification"]
#tally <- aggregate(Freq~Classification, final_small, sum)

#check whether depth correlates with total TIPs
# totalTIPs <- aggregate(Freq ~ Sample, final_data, sum)
# depth <- read.csv('depth.txt', header=F, stringsAsFactors = F)
# depth$Sample <- samplekey[match(depth$V1, samplekey$V1),]$V2
# depth$milreads <- as.numeric(gsub(".*?([0-9]+).*", "\\1", depth$V2))
# depth$totalTIPs <- totalTIPs[match(depth$Sample, totalTIPs$Sample),]$Freq
# ggplot(depth, aes(x=milreads, y=totalTIPs)) + geom_point()



# ttest <- function(famname) {
#   mini <- subset(mobile_fams, TE==famname & (Population=="S-plastotype" | Population=="D-plastotype"))
#   return(t.test(
#     subset(mini, Population=="D-plastotype")$Freq,
#     subset(mini, Population=="S-plastotype")$Freq
#   )$p.value)
# }
# mobile_fams <- subset(final_data, TE %in% subset(famtotals, Freq >= 10)$TE)
# pvals <- sapply(unique(as.character(mobile_fams$TE)), function(n) ttest(n))
# p <- p.adjust(pvals, method="fdr")
# p[which(p < 0.05)]

#Only one significant result (with any method), looks very underwhelming subset(mobile_fams, TE=="Copia-18_BD-I")



