library(ggplot2)
library(tidyr)
library(stringr)

#manually copying output from horseraceUpdatedNEW.py,
#see my notebook Feb. 2 2022
horserace <- data.frame(
  'tissue' = c('Root', 'Leaf', 'Floret', 'Callus'),
  'favorsD' = c(6488, 5898, 6704, 5691),
  'total_pairs' = c(12937, 12010, 13376, 11414)
)
#OLD RNA-SEQ DATA FROM ABR113, REVISITING FOR DISSERTATION:
#see notebook, 7 Feb. 2022
horserace <- data.frame(
  'tissue' = c('Leaf', 'Leaf', 'Leaf', 'Spike', 'Spike'),
  'favorsD' = c(6065, 6169, 6159, 9207, 8132),
  'total_pairs' = c(12247, 12272, 12263, 16384, 16133)
)
horserace$favorsS <- horserace$total_pairs - horserace$favorsD
horserace <- subset(horserace, select = -c(total_pairs) )
names(horserace) <- c('tissue', 'D', 'S')
horserace <- gather(horserace, favored_subg, num_genes, D:S, factor_key=TRUE)

########################
#FOR ABR113:
horserace$mean <- c(
  rep(mean(subset(horserace, tissue=="Leaf" & favored_subg=="D")$num_genes), times=3),
  rep(mean(subset(horserace, tissue=="Spike" & favored_subg=="D")$num_genes), times=2),
  rep(mean(subset(horserace, tissue=="Leaf" & favored_subg=="S")$num_genes), times=3),
  rep(mean(subset(horserace, tissue=="Spike" & favored_subg=="S")$num_genes), times=2)
  )
horserace$sd <- c(
  rep(sd(subset(horserace, tissue=="Leaf" & favored_subg=="D")$num_genes), times=3),
  rep(sd(subset(horserace, tissue=="Spike" & favored_subg=="D")$num_genes), times=2),
  rep(sd(subset(horserace, tissue=="Leaf" & favored_subg=="S")$num_genes), times=3),
  rep(sd(subset(horserace, tissue=="Spike" & favored_subg=="S")$num_genes), times=2)
)
horserace <- subset(horserace, select = -c(num_genes) )
horserace <- horserace[!duplicated(horserace),]
########################

#tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb2_horserace.tiff", units="in", width=6.5, height=3, res=300)
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb9_ABR113_horserace.tiff", units="in", width=5, height=3, res=300)
#Bhyb26:
#ggplot(horserace, aes(x=tissue, y=num_genes, fill=favored_subg))+
#ABR113:
ggplot(horserace, aes(x=tissue, y=mean, fill=favored_subg))+
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


#manually copying output from analyze_counts_by_subg.py,
#see my notebook Nov.29 2021
countdat <- data.frame(
  'tissue' = c('Root', 'Leaf', 'Floret', 'Callus'),
  'D' = c(502979, 622342, 709650, 221374),
  'S' = c(485971, 571728, 819356, 205699)
)
countdat <- gather(countdat, subg, num_counts, D:S, factor_key=TRUE)

##################
#OLD RNA-SEQ DATA FROM ABR113, REVISITING FOR DISSERTATION:
countdat <- data.frame(
  'tissue' = c('Leaf', 'Leaf', 'Leaf', 'Spike', 'Spike'),
  'D' = c(15165822, 18195689, 17815259, 13150381, 12309994),
  'S' = c(15473884, 18337787, 18075138, 12816582, 11339499)
)
countdat <- gather(countdat, subg, num_counts, D:S, factor_key=TRUE)
countdat$sample <- c('Leaf 1', 'Leaf 2', 'Leaf 3', 'Spike 1', 'Spike 2', 'Leaf 1', 'Leaf 2', 'Leaf 3', 'Spike 1', 'Spike 2')
#ABR113:
# wc -l BhDprimTrlengths.txt 
# 37711 BhDprimTrlengths.txt
# wc -l BhSprimTrlengths.txt 
# 32449 BhSprimTrlengths.txt
# 32449/(32449+37711) = 0.4625
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb9_subgcount_ABR113_leaf.tiff", units="in", width=6, height=4, res=300)
ggplot(subset(countdat, tissue=="Leaf"), aes(x=sample, y=num_counts, fill=subg))+
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_hline(yintercept = 0.4625, color='red', size=1)+ 
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
tiff("/home/virginia/Documents/school/vogelLab/notebook/2022/feb9_subgcount_ABR113_spike.tiff", units="in", width=5, height=4, res=300)
ggplot(subset(countdat, tissue=="Spike"), aes(x=sample, y=num_counts, fill=subg))+
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_hline(yintercept = 0.4625, color='red', size=1)+ 
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
###############





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


#Adding some code Feb. 2 2022: So, homeologs in general are not biased toward one subgenome.
#What about the "significantly" biased homeolog pairs? I put "significantly" in quotes because
#this is not a proper DEG analysis. See my notebook Jan 21 2022.

#Lifting some code from TPMs.R.
#Not  including pseudogenes in this analysis.
TPMs <- read.table('/home/virginia/Documents/school/vogelLab/notebook/2022/all_TPMs_w_pseudogenes.txt', sep='\t', header = TRUE, stringsAsFactors = FALSE)
TPMs$Genetype <- sapply(TPMs$Gene, function(x) ifelse(startsWith(x, "BrBhy"), "Conserved genes", "Putative pseudogenes"))
homeo <- read.table(file="/home/virginia/Documents/school/vogelLab/notebook/2021/Bhyb26_homeologs.tsv", sep = "\t", header = T, stringsAsFactors = F)
names(homeo) <- c("D", "S")

#log(b/a), which is the same as log(b) - log(a), where a and b are the homeolog expression levels:
#Smith et al. 2019 refer to this as "a dimensionless quantity with HEB=0 indicating no bias. 
#If one uses the base 2 logarithm, HEB=−3 indicates 8-fold bias towards homeolog A."
#In this case, we are doing log2(d/s), so 0 means no bias, negative means biased toward S,
#and positive means biased toward D

D_TPM_root <- sapply(homeo$D, function(x) TPMs[match(x, TPMs$Gene),]$Root)
S_TPM_root <- sapply(homeo$S, function(x) TPMs[match(x, TPMs$Gene),]$Root)
HEB_root <- mapply(function(d, s) {log2((d+0.1)/(s+0.1))}, D_TPM_root, S_TPM_root)
HEB_root_df <- data.frame(value=unlist(unname(HEB_root)), stringsAsFactors = F)
#root_nozero <- HEB_root[HEB_root!=0]
#root_nozero_df <- data.frame(value=unlist(unname(root_nozero)), stringsAsFactors = F)
#ggplot(root_nozero_df, aes(x=value)) + geom_histogram(binwidth=0.5) 

D_TPM_leaf <- sapply(homeo$D, function(x) TPMs[match(x, TPMs$Gene),]$Leaf)
S_TPM_leaf <- sapply(homeo$S, function(x) TPMs[match(x, TPMs$Gene),]$Leaf)
HEB_leaf <- mapply(function(d, s) {log2((d+0.1)/(s+0.1))}, D_TPM_leaf, S_TPM_leaf)
HEB_leaf_df <- data.frame(value=unlist(unname(HEB_leaf)), stringsAsFactors = F)
#leaf_nozero <- HEB_leaf[HEB_leaf!=0]
#leaf_nozero_df <- data.frame(value=unlist(unname(leaf_nozero)), stringsAsFactors = F)
#ggplot(leaf_nozero_df, aes(x=value)) + geom_histogram(binwidth=0.5) 

D_TPM_floret <- sapply(homeo$D, function(x) TPMs[match(x, TPMs$Gene),]$Floret)
S_TPM_floret <- sapply(homeo$S, function(x) TPMs[match(x, TPMs$Gene),]$Floret)
HEB_floret <- mapply(function(d, s) {log2((d+0.1)/(s+0.1))}, D_TPM_floret, S_TPM_floret)
HEB_floret_df <- data.frame(value=unlist(unname(HEB_floret)), stringsAsFactors = F)

D_TPM_callus <- sapply(homeo$D, function(x) TPMs[match(x, TPMs$Gene),]$Callus)
S_TPM_callus <- sapply(homeo$S, function(x) TPMs[match(x, TPMs$Gene),]$Callus)
HEB_callus <- mapply(function(d, s) {log2((d+0.1)/(s+0.1))}, D_TPM_callus, S_TPM_callus)
HEB_callus_df <- data.frame(value=unlist(unname(HEB_callus)), stringsAsFactors = F)

homeoBias <- cbind(homeo, HEB_root_df, HEB_leaf_df, HEB_floret_df, HEB_callus_df)

get_lower_tail <- function(myvec) {
  return(quantile(myvec, probs=seq(0, 1, 1/200), na.rm=TRUE)[["2.5%"]])
}
get_upper_tail <- function(myvec) {
  return(quantile(myvec, probs=seq(0, 1, 1/200), na.rm=TRUE)[["97.5%"]])
}
get_sig <- function(colnum) {
  return(
    c(which(homeoBias[,colnum] < get_lower_tail(homeoBias[,colnum])),
      which(homeoBias[,colnum] > get_upper_tail(homeoBias[,colnum])))
  )
}

#These are just indices of homeoBias
sig_genes_root <- get_sig(3)
sig_genes_leaf <- get_sig(4)
sig_genes_floret <- get_sig(5)
sig_genes_callus <- get_sig(6)

sig_in_any_expt <- unique(c(
  sig_genes_root,
  sig_genes_leaf,
  sig_genes_floret,
  sig_genes_callus
))

sig_in_all_expts <- Reduce(intersect, list(
  sig_genes_root,
  sig_genes_leaf,
  sig_genes_floret,
  sig_genes_callus
))

conditionally_biased_pairs <- homeoBias[sig_in_any_expt,]
consistently_biased_pairs <- homeoBias[sig_in_all_expts,]
nrow(conditionally_biased_pairs)
nrow(consistently_biased_pairs)
# > nrow(conditionally_biased_pairs)
# [1] 2460
# > nrow(consistently_biased_pairs)
# [1] 212

homeoBias$root_favored <- NA
homeoBias$leaf_favored <- NA
homeoBias$floret_favored <- NA
homeoBias$callus_favored <- NA

homeoBias[sig_genes_root,]$root_favored <- ifelse(homeoBias[sig_genes_root,3] < 0, "S", "D")
homeoBias[sig_genes_leaf,]$leaf_favored <- ifelse(homeoBias[sig_genes_leaf,4] < 0, "S", "D")
homeoBias[sig_genes_floret,]$floret_favored <- ifelse(homeoBias[sig_genes_floret,5] < 0, "S", "D")
homeoBias[sig_genes_callus,]$callus_favored <- ifelse(homeoBias[sig_genes_callus,6] < 0, "S", "D")

#Interestingly, 100% of homeolog pairs that were biased in all experiments were biased
#in the same direction in all experiments. (though n=212 so small sample size)
homeoBias[rownames(consistently_biased_pairs),c(7, 8, 9, 10)]

#To do: what proportion of consistently biased pairs favored D? S?
nrow(
  subset(homeoBias[rownames(consistently_biased_pairs),c(7, 8, 9, 10)], 
       root_favored=="S" & leaf_favored =="S" & floret_favored=="S" & callus_favored=="S")
)
nrow(
  subset(homeoBias[rownames(consistently_biased_pairs),c(7, 8, 9, 10)], 
         root_favored=="D" & leaf_favored =="D" & floret_favored=="D" & callus_favored=="D")
)
#95 pairs favored S, 117 pairs favored D
binom.test(117, 212, p=0.5, alternative = c("two.sided"))
# Exact binomial test
# 
# data:  117 and 212
# number of successes = 117, number of trials = 212, p-value = 0.149
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#   0.4822611 0.6200411
# sample estimates:
#   probability of success 
# 0.5518868 

#To do: how many cases of significant HEB favored D? S?
nrow(subset(homeoBias, root_favored=="S"))
nrow(subset(homeoBias, root_favored=="D"))
nrow(subset(homeoBias, leaf_favored=="S"))
nrow(subset(homeoBias, leaf_favored=="D"))
nrow(subset(homeoBias, floret_favored=="S"))
nrow(subset(homeoBias, floret_favored=="D"))
nrow(subset(homeoBias, callus_favored=="S"))
nrow(subset(homeoBias, callus_favored=="D"))

#BAH! Of course it's the same number in all of them (502). Because we are taking the top and bottom
#2.5%, we are taking the top n S-favoring pairs and the top n D-favoring pairs.
# > length(sig_genes_root)
# [1] 1004
# > length(sig_genes_leaf)
# [1] 1004
# > nrow(homeo)
# [1] 20066
  # > .025*20066
# [1] 501.65

#I guess another way to look at this is to ask, what is the threshold for significance for each subgenome
#in each tissue?

get_lower_cutoff <- function(myvec) {
  return(quantile(myvec, probs=seq(0, 1, 1/200), na.rm=TRUE)[["2.5%"]])
}
get_upper_cutoff <- function(myvec) {
  return(quantile(myvec, probs=seq(0, 1, 1/200), na.rm=TRUE)[["97.5%"]])
}

#root S cutoff
get_lower_cutoff(homeoBias[,3])
#root D cutoff
get_upper_cutoff(homeoBias[,3])
#leaf S cutoff
get_lower_cutoff(homeoBias[,4])
#leaf D cutoff
get_upper_cutoff(homeoBias[,4])
#floret S cutoff
get_lower_cutoff(homeoBias[,5])
#floret D cutoff
get_upper_cutoff(homeoBias[,5])
#callus S cutoff
get_lower_cutoff(homeoBias[,6])
#callus D cutoff
get_upper_cutoff(homeoBias[,6])

# > #root S cutoff
#   > get_lower_cutoff(homeoBias[,3])
# [1] -6.045022
# > #root D cutoff
#   > get_upper_cutoff(homeoBias[,3])
# [1] 6.093473
# > #leaf S cutoff
#   > get_lower_cutoff(homeoBias[,4])
# [1] -6.108749
# > #leaf D cutoff
#   > get_upper_cutoff(homeoBias[,4])
# [1] 6.057693
# > #floret S cutoff
#   > get_lower_cutoff(homeoBias[,5])
# [1] -5.60152
# > #floret D cutoff
#   > get_upper_cutoff(homeoBias[,5])
# [1] 5.735106
# > #callus S cutoff
#   > get_lower_cutoff(homeoBias[,6])
# [1] -7.005793
# > #callus D cutoff
#   > get_upper_cutoff(homeoBias[,6])
# [1] 7.081963



  