library(ggplot2)
library(RColorBrewer)
options(scipen=999)

get_sum <- function(submat){
  return(sapply(seq(nrow(submat)), function(n) sum(as.numeric(submat[n,]))))
}

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
test <- c("gray70", "#C269CF", "#1FA9A3", "#FEE227", "#0F61AF", "#FE036A", "#20B318", "#613E27")
test2 <- c("gray70", colorRampPalette(c("tomato3", "gold", "seagreen3", "deepskyblue4","purple3", "plum2"))(7))


setwd('/home/virginia/Documents/school/vogelLab/notebook/2022')
ABR113D_subg_size <- 269318743
ABR113S_subg_size <- 239897624
Bhyb26D_subg_size <-277878319
Bhyb26S_subg_size <- 249002919
dis_genome_size <- 271067295
sta_genome_size <- 233737016

ABR113_BhD_dat <- read.csv('bp_bychrom.BhD.txt_ABR113', sep='\t', header=TRUE, stringsAsFactors=FALSE) #file produced by getTEstats.py in step6 of my pipeline
ABR113_BhD_dat$sum <- get_sum(ABR113_BhD_dat[,which(startsWith(names(ABR113_BhD_dat), "BhD"))])
ABR113_BhS_dat <- read.csv('bp_bychrom.BhS.txt_ABR113', sep='\t', header=TRUE, stringsAsFactors=FALSE) #file produced by getTEstats.py in step6 of my pipeline
ABR113_BhS_dat$sum <- get_sum(ABR113_BhS_dat[,which(startsWith(names(ABR113_BhS_dat), "BhS"))])

Bhyb26_BhD_dat <- read.csv('bp_bychrom.BhD.txt_Bhyb26', sep='\t', header=TRUE, stringsAsFactors=FALSE) #file produced by getTEstats.py in step6 of my pipeline
Bhyb26_BhD_dat$sum <- get_sum(Bhyb26_BhD_dat[,which(startsWith(names(Bhyb26_BhD_dat), "BhD"))])
Bhyb26_BhS_dat <- read.csv('bp_bychrom.BhS.txt_Bhyb26', sep='\t', header=TRUE, stringsAsFactors=FALSE) #file produced by getTEstats.py in step6 of my pipeline
Bhyb26_BhS_dat$sum <- get_sum(Bhyb26_BhS_dat[,which(startsWith(names(Bhyb26_BhS_dat), "BhS"))])

#BTW I split those files^ by subgenome like so:
# setwd('/home/virginia/Documents/school/vogelLab/notebook/2022')
# dat <- read.csv('bp_bychrom.txt_Bhyb26', sep='\t', header=TRUE, stringsAsFactors=FALSE) #file produced by getTEstats.py in step6 of my pipeline
# Ddat <- dat[,c(1:6)]
# Sdat <- dat[,c(1, 7:16)]
# write.table(Ddat, file='bp_bychrom.BhD.txt_Bhyb26', row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
# write.table(Sdat, file='bp_bychrom.BhS.txt_Bhyb26', row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")

dis_dat <- read.csv('bp_bychrom.txt_Bd21', sep='\t', header=TRUE, stringsAsFactors=FALSE)
sta_dat <- read.csv('bp_bychrom.txt_ABR114', sep='\t', header=TRUE, stringsAsFactors=FALSE)

dis_dat$sum <- sapply(seq(nrow(dis_dat)), function(n) sum(as.numeric(dis_dat[n,2:ncol(dis_dat)])))
sta_dat$sum <- sapply(seq(nrow(sta_dat)), function(n) sum(as.numeric(sta_dat[n,2:ncol(sta_dat)])))

TEsums <- data.frame(Genome = c(rep("ABR113D", times=nrow(ABR113_BhD_dat)), rep("ABR113S", times=nrow(ABR113_BhS_dat)), rep("Bhyb26D", times=nrow(Bhyb26_BhD_dat)), rep("Bhyb26S", times=nrow(Bhyb26_BhS_dat)), rep("Bd21", times=nrow(dis_dat)), rep("ABR114", times=nrow(sta_dat))), 
                        Classification = c(ABR113_BhD_dat$Classification, ABR113_BhS_dat$Classification, Bhyb26_BhD_dat$Classification, Bhyb26_BhS_dat$Classification, dis_dat$Classification, sta_dat$Classification),
                        Value = c(ABR113_BhD_dat$sum, ABR113_BhS_dat$sum, Bhyb26_BhD_dat$sum, Bhyb26_BhS_dat$sum, dis_dat$sum, sta_dat$sum))

nonTEsums <- data.frame(Genome = c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bd21", "ABR114"), 
                    Classification = "non-TE genome space", 
                    Value = c(ABR113D_subg_size-sum(ABR113_BhD_dat$sum), ABR113S_subg_size-sum(ABR113_BhS_dat$sum), Bhyb26D_subg_size-sum(Bhyb26_BhD_dat$sum), Bhyb26S_subg_size-sum(Bhyb26_BhS_dat$sum), dis_genome_size-sum(dis_dat$sum), sta_genome_size-sum(sta_dat$sum)))

#Putting all the MITEs into one category since it's visually distracting to keep them separate
mites <- subset(TEsums, grepl( "MITE", Classification, fixed = TRUE))
newmites <- data.frame(Genome = c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bd21", "ABR114"), 
                       Classification = "DNA_transposon_MITE",
                       Value = c(sum(subset(mites, Genome=="ABR113D")$Value), 
                                 sum(subset(mites, Genome=="ABR113S")$Value),
                                 sum(subset(mites, Genome=="Bhyb26D")$Value), 
                                 sum(subset(mites, Genome=="Bhyb26S")$Value),
                                 sum(subset(mites, Genome=="Bd21")$Value),
                                 sum(subset(mites, Genome=="ABR114")$Value)))

clean <- subset(TEsums, !(grepl( "MITE", Classification, fixed = TRUE)) & !(grepl( "could_not_classify", Classification, fixed = TRUE)))
all_dat <- rbind(nonTEsums, clean, newmites)
all_dat$Classification <- factor(all_dat$Classification, levels=c(
  "non-TE genome space",
  unique(ABR113_BhD_dat$Classification)[which(!(grepl( "MITE", unique(ABR113_BhD_dat$Classification), fixed = TRUE)) & !(grepl( "could_not_classify", unique(ABR113_BhD_dat$Classification), fixed = TRUE)))],
  "DNA_transposon_MITE"))
all_dat$Genome <- factor(as.character(all_dat$Genome), levels=c("Bd21", "ABR113D", "Bhyb26D", "Bhyb26S", "ABR113S", "ABR114"))

#In fact, let's put all DNA TEs in one category except the two biggest, MITE and CACTA
DNATEs <- subset(all_dat,
  grepl("DNA_transposon", Classification, fixed = TRUE) & 
    !(grepl("CACTA", Classification, fixed=T)) & 
    !(grepl("MITE", Classification, fixed=T)))

newDNA <- data.frame(Genome = c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bd21", "ABR114"), 
                       Classification = "Other_DNA_TE",
                       Value = c(sum(subset(DNATEs, Genome=="ABR113D")$Value), 
                                 sum(subset(DNATEs, Genome=="ABR113S")$Value),
                                 sum(subset(DNATEs, Genome=="Bhyb26D")$Value), 
                                 sum(subset(DNATEs, Genome=="Bhyb26S")$Value),
                                 sum(subset(DNATEs, Genome=="Bd21")$Value),
                                 sum(subset(DNATEs, Genome=="ABR114")$Value)))

noDNA <- subset(all_dat, !grepl("DNA_transposon", Classification, fixed = TRUE))
bigDNA <- subset(all_dat, (grepl("CACTA", Classification, fixed=T) | grepl("MITE", Classification, fixed=T)) )
toPlot <- rbind(noDNA, bigDNA, newDNA)

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/mar28_subg_overview.tiff", units="in", width=4.5, height=4.25, res=350)
ggplot(data=toPlot, aes(y=Value/1000000, x=Genome, fill=Classification)) +
  geom_bar(position="stack", stat="identity" )+
  ylab("Genome Space (Mb)")+
  scale_y_continuous(limits=c(0, 300), breaks=c(50,100,150,200,250), expand = c(0, 0))+
  theme(
    legend.position = "right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    #legend.spacing.y = unit(0.5, 'cm'),
    #legend.key.height = unit(0.5, 'cm'),
    #legend.key.width = unit(0.5, 'cm'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=15, angle=-60, hjust=0.05, vjust=0.5),
    axis.text.y = element_text(size=15),
    panel.border=element_rect(colour="gray", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
    panel.grid.minor = element_line(size = 0.2, linetype = 'solid', color = "gray90"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  scale_fill_manual(values = test, labels=c("non-TE", "LTR-RT RLC", "LTR-RT RLG", "LTR-RT RLX","non-LTR RT","CACTA", "MITE", "Other DNA TE")) 


dev.off()

# ABR113.chrom.sizes:
# BhD1	73190852
# BhD2	59603763
# BhD3	59465344
# BhD4	48382457
# BhD5	28676327
# BhS1	30611413
# BhS10	19337815
# BhS2	28347457
# BhS3	26199158
# BhS4	25450061
# BhS5	23650338
# BhS6	22645271
# BhS7	21641948
# BhS8	20999971
# BhS9	21014192

# cat Bhyb26/Bhyb26.chrom.sizes
# BhD1	78242934
# BhD2	60904068
# BhD3	62918427
# BhD4	48441464
# BhD5	27371426
# BhS1	31901534
# BhS2	30385391
# BhS3	26678100
# BhS4	26298841
# BhS5	23487869
# BhS6	23238163
# BhS7	21996878
# BhS8	22147626
# BhS9	22306114
# BhS10	20562403

# 78242934+60904068+62918427+48441464+27371426
# 31901534+30385391+26678100+26298841+23487869+23238163+21996878+22147626+22306114+20562403

# Bd21.chroms.sizes
# Bd1	75071545
# Bd2	59130575
# Bd3	59640145
# Bd4	48594894
# Bd5	28630136
#75071545+59130575+59640145+48594894+28630136

# ABR114.chrom.sizes
# Chr01   30108798
# Chr02   27797795
# Chr03   25300955
# Chr04   24648699
# Chr05   23060899
# Chr06   21822193
# Chr07   20898793
# Chr08   20713557
# Chr09   20580347
# Chr10   18804980
#30108798+27797795+25300955+24648699+23060899+21822193+20898793+20713557+20580347+18804980



#clean2 <- rbind(clean[c(1:12),],newmites[1,],clean[13:nrow(clean),])
#TEsums_censored <- rbind(clean2[-nrow(clean2),],newmites[2,],clean2[nrow(clean2),])
#TEsums_censored$Classification <- as.character(TEsums_censored$Classification)
#final <- rbind(nonTEsums, TEsums_censored)

# justTEs <- subset(all_dat, !(grepl( "non-TE genome space", Classification, fixed = TRUE)))
# ggplot(data=justTEs, aes(y=Value/1000000, x=Genome, fill=Classification)) +
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_manual(values = mycolors) +
#   ylab("Genome Space (Mb)")+
#   scale_y_continuous(limits=c(0, 280), breaks=c(50,100,150,200,250))+
#   theme(legend.title=element_blank())


#   
# ggplot(justTEs, aes(fill=Classification, y=Value, x=Genome)) + 
#   geom_bar(position="fill", stat="identity")+
#   scale_fill_manual(values = mycolors) +
#   ylab("Proportion of TE space")+
#   theme(legend.title=element_blank())
