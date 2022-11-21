library(ggplot2)
options(scipen=999)

# "digital chromosome painting"
# This code is super messy because it took me SO LONG to get this plot right
# Takes a csv file with columns: Chr,Start,End,Value
# As you can see from the comments, I have used this script
# to plot many types of data over the years

#setwd('/home/virginia/Desktop')
setwd('/home/virginia/Documents/school/vogelLab/notebook/2021')
chromsizesfile <- 'ABR113.chrom.sizes'
#chromsizesfile <- 'Bhyb26.chrom.sizes'
#chromsizesfile <- 'Bd21.chrom.sizes'
#chromsizesfile <- 'ABR114.chrom.sizes'
#densityfile <- 'rnd-4_family-1045_ABR113_density_by_bin.tsv'
#densityfile <- 'TEdensity.Bhyb26.200kb.csv'
#densityfile <- 'pacbio_BAC_CSNPB_frag_densities.csv'
#densityfile <- 'pacbio_BAC_CSNPB_frag_densities_BHYB26.csv'
#densityfile <- 'subgenome_specific_TEs.csv_ABR113' # see scripts/densities/getFragDensity.R and get_subg_specific_TE_density.R
#densityfile <- 'subgenome_specific_TEs.Bhyb26.80perc.csv'
#densityfile <- 'lost_gene_density.csv' #from my script pseudogene_density.R CURRENTLY IN MANUSCRIPT 
#densityfile <- 'TEs_common_to_Bd21Bd11Bd301.csv'
#densityfile <- 'TEs_private_to_Bd21_vs_ABR114.csv'
#densityfile <- 'TEs_common_to_Bd21ABR114.csv'
#densityfile <- 'TEs_private_to_Bd21.csv'
#densityfile <- 'TEs_private_to_ABR114_vs_Bd21.csv'
#densityfile <- 'TEs_common_to_ABR114Bd21.csv'
#densityfile <- 'TEs_Bd21.csv'
#densityfile <- 'TEs_ABR113.csv'
densityfile <- 'pacbioBACs/CSNPB_1_CSNPC_3.bins'

dat <- read.csv(densityfile, header=TRUE, stringsAsFactors = FALSE)
windowSize <- dat[1,'End']
dat <- dat[which(dat$End-dat$Start==windowSize-1),] #remove the last bin on each chromosome, which is smaller than the windowsize
#^actually we are saying, "keep only those bins that are windowSize-1 bp long"

chrdat <- read.csv(chromsizesfile, sep='\t', header=FALSE, stringsAsFactors = FALSE)
chrs_ordered <- c("BhS10","BhS9","BhS8","BhS7","BhS6","BhS5","BhS4","BhS3","BhS2","BhS1","BhD5","BhD4","BhD3","BhD2","BhD1")
#chrs_ordered <- c("Chr10","Chr09","Chr08","Chr07","Chr06","Chr05","Chr04","Chr03","Chr02","Chr01")
#chrs_ordered <- c("Bd1", "Bd2", "Bd3", "Bd4", "Bd5")
dat$bin_num <- unlist(sapply(chrdat$V1, function(x) seq(unname(nrow(subset(dat, Chr==x))))))

increment <- 1000000/windowSize
#create a df with the dimensions of each chromosome
df_box <- data.frame(Chr=chrs_ordered)
df_box$x_min <- rep(0, times=length(chrs_ordered))
df_box$x_max <- as.numeric(unlist(sapply(chrs_ordered, function(x) max(subset(dat, Chr==x)$bin_num+1)/increment)))
df_box$y_min <- sapply(seq(0,length(chrs_ordered)-1), function(x) x+0.5)
df_box$y_max <- sapply(seq(1,length(chrs_ordered)), function(x) x+0.5)

#adding some space between the boxes
to_add <- sapply(seq(nrow(df_box)), function(x) x*0.5-0.5)
df_box$y_min <- df_box$y_min + to_add
df_box$y_max <- df_box$y_max + to_add
df_box$y_mean <-sapply(seq(nrow(df_box)), function(x) mean(c(df_box[x,'y_min'], df_box[x,'y_max'])) )

dat$y_min <- df_box[match(dat$Chr, df_box$Chr),'y_min']
dat$y_max <- df_box[match(dat$Chr, df_box$Chr),'y_max']
dat$y_mean <- df_box[match(dat$Chr, df_box$Chr),'y_mean']
dat$x_min <- rep(0, times=nrow(dat))
dat$x_max <- df_box[match(dat$Chr, df_box$Chr),'x_min']

tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/apr6_pacbioBACs_CSNPBbCSNPCg.tiff", units="in", width=4, height=5, res=400)
ggplot(dat, aes(x=bin_num/increment, y=y_mean, width=1, height=1))+
#geom_tile(aes(fill = Pseudogenes*100))+
geom_tile(aes(fill = factor(Value)))+
  theme(
    legend.position="none",
    #legend.title=element_text(size=14), 
    #legend.text=element_text(size=14),
    #legend.margin=margin(0,0,0,0),
    #legend.box.margin=margin(-5,-5,-5,-5),
    panel.background = element_rect(fill = "white", colour = "grey50"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
    aspect.ratio = 1)+
  scale_x_continuous (name = "Location (Mbp)") +
  scale_y_continuous(name = "Chromosome", breaks=df_box$y_mean, labels=chrs_ordered)+
  #scale_fill_manual(values=c("white", "blue", "red", "green"), breaks=c(0, 1, 2, 3), labels=c(0, 1, 2, 3))+
  scale_fill_gradient2(low = "white", high = "black")+
  #scale_fill_gradient2(low = "white", high = "red")+ #guide =  guide_legend("TE Count")) + #guide_legend("TE Density (% of bin)")) +
  geom_rect(data = df_box[,c(2,3,4,5)], 
            aes(xmin = x_min, xmax = x_max,
                ymin = y_min, ymax = y_max),
            inherit.aes = FALSE,
            colour = "black", fill=NA, size = 0.5)
dev.off()


