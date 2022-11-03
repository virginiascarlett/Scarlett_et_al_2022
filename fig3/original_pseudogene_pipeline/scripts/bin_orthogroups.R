library(ggplot2)
library(data.table)
library(UpSetR)

get_og_comp <- function(myog) {
  return(DT[synOg==myog]$genome)
}
#has_syn_orth checks:
#if the gene is from ABR113D, is at least one ABR113S gene in this syntenic orthogroup?
#if the gene is from Bhyb26S, is at least one Bhyb26D gene in this syntenic orthogroup?
#etc.
#if the gene is from Bdistachyon, is at least one Bstacei gene in this syntenic orthogroup?
#if the gene is from Bstacei, is at least one Bdistachyon gene in this syntenic orthogroup?
has_syn_orth <- function(geneID){
  generow <- DT[id==geneID,]
  genome_orig <- generow$'genome'
  ogcomp <- ogkey[[generow$'synOg']] #returns a vector that looks like e.g. "ABR113D"     "ABR113S"     "Bdistachyon" "Bhyb26S"     "Bstacei"     "Osativa"
  if (startsWith(genome_orig, 'ABR113')) {
    if (genome_orig=="ABR113D") {return(ifelse(length(which(ogcomp=="ABR113S")>0 ), TRUE, FALSE))}
  } else {
    return(ifelse(length(which(ogcomp=="ABR113D")>0 ), TRUE, FALSE))
  }
  if (startsWith(genome_orig, 'Bhyb26')) {
    if (genome_orig=='Bhyb26D') {return(ifelse(length(which(ogcomp=="Bhyb26S")>0 ), TRUE, FALSE))}
  } else {
    return(ifelse(length(which(ogcomp=="Bhyb26D")>0 ), TRUE, FALSE))
  }
  if (startsWith(genome_orig, 'Bd21')) {
    return(ifelse(length(which(ogcomp=="Bstacei")>0 ), TRUE, FALSE))
  }
  if (startsWith(genome_orig, 'Bstacei')) {
    return(ifelse(length(which(ogcomp=="Bd21")>0 ), TRUE, FALSE))
  }
  if (!any(grepl(genome_orig, ogcomp))) { #for some reason this line isn't working, doesn't print anything if given a rice gene
    message("ERROR: geneID is not from Bd21, Bstacei, ABR113 or Bhyb26")
  }
}

setwd("/home/virginia/Documents/school/vogelLab/notebook/2021")
dat <- read.csv('/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes/results/gffWithOgs.txt', sep='\t', stringsAsFactors = FALSE, header = TRUE)
DT <- as.data.table(dat)     


#takes about a minute
ogkey <- lapply(unique(dat$synOg), function(myog) get_og_comp(myog))
names(ogkey) <- unique(dat$synOg)

#Check what proportion of genes has a syntenic ortholog
#takes about a minute EACH
#skip if you don't need this
#Expects Bdistachyon genome to be "Bd21"
# syn_orth_Bhyb26D <- sapply(subset(dat, genome=='Bhyb26D')$id, function(mygene) has_syn_orth(mygene) )
# syn_orth_Bhyb26S <- sapply(subset(dat, genome=='Bhyb26S')$id, function(mygene) has_syn_orth(mygene) )
# syn_orth_ABR113D <- sapply(subset(dat, genome=='ABR113D')$id, function(mygene) has_syn_orth(mygene) )
# syn_orth_ABR113S <- sapply(subset(dat, genome=='ABR113S')$id, function(mygene) has_syn_orth(mygene) )
# syn_orth_Bd21 <- sapply(subset(dat, genome=='Bdistachyon')$id, function(mygene) has_syn_orth(mygene) )
# syn_orth_Bstacei <- sapply(subset(dat, genome=='Bstacei')$id, function(mygene) has_syn_orth(mygene) )
# 
# sum(syn_orth_Bhyb26D)/length(syn_orth_Bhyb26D)
# sum(syn_orth_Bhyb26S)/length(syn_orth_Bhyb26S)
# sum(syn_orth_ABR113D)/length(syn_orth_ABR113D)
# sum(syn_orth_ABR113S)/length(syn_orth_ABR113S)
# sum(syn_orth_Bd21)/length(syn_orth_Bd21)
# sum(syn_orth_Bstacei)/length(syn_orth_Bstacei)


#Check which syntenic orthogroups have this combination of genomes represented
#Does not depend on the earlier syn_orth_ stuff so you can just skip to here 
conserved <- lapply(ogkey, function(genomes) all(c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bstacei", "Bdistachyon", "Osativa") %in% genomes))
abs_Bhyb26D <- lapply(ogkey, function(genomes) (all(c("ABR113D", "ABR113S", "Bhyb26S", "Bstacei", "Bdistachyon", "Osativa") %in% genomes) & !("Bhyb26D" %in% genomes)))
abs_Bhyb26S <- lapply(ogkey, function(genomes) (all(c("ABR113D", "ABR113S", "Bhyb26D", "Bstacei", "Bdistachyon", "Osativa") %in% genomes) & !("Bhyb26S" %in% genomes)))
abs_ABR113D <- lapply(ogkey, function(genomes) (all(c("ABR113S", "Bhyb26D", "Bhyb26S", "Bstacei", "Bdistachyon", "Osativa") %in% genomes) & !("ABR113D" %in% genomes)))
abs_ABR113S <- lapply(ogkey, function(genomes) (all(c("ABR113D", "Bhyb26D", "Bhyb26S", "Bstacei", "Bdistachyon", "Osativa") %in% genomes) & !("ABR113S" %in% genomes)))
abs_Bd21 <- lapply(ogkey, function(genomes) (all(c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bstacei", "Osativa") %in% genomes) & !("Bdistachyon" %in% genomes)))
abs_Bstacei <- lapply(ogkey, function(genomes) (all(c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bdistachyon", "Osativa") %in% genomes) & !("Bstacei" %in% genomes)))
abs_rice <- lapply(ogkey, function(genomes) (all(c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bdistachyon", "Bstacei") %in% genomes) & !("Osativa" %in% genomes)))

input = c(
  'Bhyb26D&Bhyb26S&ABR113D&ABR113S&Bd21&ABR114'=sum(unlist(abs_rice)), 
  'Osativa&Bhyb26S&ABR113D&ABR113S&Bd21&ABR114'=sum(unlist(abs_Bhyb26D)),
  'Osativa&Bhyb26D&ABR113D&ABR113S&Bd21&ABR114'=sum(unlist(abs_Bhyb26S)),
  'Osativa&Bhyb26D&Bhyb26S&ABR113S&Bd21&ABR114'=sum(unlist(abs_ABR113D)),
  'Osativa&Bhyb26D&Bhyb26S&ABR113D&Bd21&ABR114'=sum(unlist(abs_ABR113S)),
  'Osativa&Bhyb26D&Bhyb26S&ABR113D&ABR113S&ABR114'=sum(unlist(abs_Bd21)),
  'Osativa&Bhyb26D&Bhyb26S&ABR113D&ABR113S&Bd21'=sum(unlist(abs_Bstacei))
)

tiff("nov2_upset.tiff", units="in", width=7, height=5.25, res=300)
upset(fromExpression(input),
      nsets=7,
      mb.ratio = c(0.60, 0.40), 
      order.by = "freq",
      decreasing=T,
      point.size = 2.8,
      line.size = 0.8,
      text.scale=c(2.3, 2.25, 0.5, 0.5, 2.25, 3.25),
      #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      #text.scale=2.25,
      mainbar.y.label = "Number of Orthogroups",
      set_size.show=F,
      set_size.scale_max=0.5)
dev.off()

############## Get homeologs for the two polyploids ##############  
has_homeos <- function(string1, string2, testvec) {
  return(sum(testvec %in% string1)== 1 & sum(testvec %in% string2) == 1) #true or false: one and only one occurrence of each string
}
get_homeos <- function(string1, string2, mysynOg) {
  minidf <- DT[which(DT$synOg==mysynOg),]
  homeo1 <- minidf[which(minidf$genome==string1),]$id
  homeo2 <- minidf[which(minidf$genome==string2),]$id
  return(c(homeo1, homeo2))
}

synOGs_w_ABR113_homeos <- unname(which(unlist(lapply(ogkey, function(genomes) has_homeos("ABR113D", "ABR113S", genomes)))))
synOGs_w_Bhyb26_homeos <- unname(which(unlist(lapply(ogkey, function(genomes) has_homeos("Bhyb26D", "Bhyb26S", genomes)))))

ABR113homeos <- t(as.data.frame(lapply(synOGs_w_ABR113_homeos, function(n) get_homeos("ABR113D", "ABR113S", n)), stringsAsFactors = F))
Bhyb26homeos <- t(as.data.frame(lapply(synOGs_w_Bhyb26_homeos, function(n) get_homeos("Bhyb26D", "Bhyb26S", n)), stringsAsFactors = F))
write.table(ABR113homeos, file="/home/virginia/Documents/school/vogelLab/notebook/2021/ABR113_homeologs.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Bhyb26homeos, file="/home/virginia/Documents/school/vogelLab/notebook/2021/Bhyb26_homeologs.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
##################################################################




############## Test for significance of biased gene loss ############## 
get_total_genes_in_subg <- function(mygenome) {
  return(
    length(unique(
      DT[which(DT$genome==mygenome),]$id
    ))
  )
}

#Bhyb26
contingency_table <- data.frame("D"=c(sum(unlist(abs_Bhyb26D)), get_total_genes_in_subg("Bhyb26D")), "S"=c(sum(unlist(abs_Bhyb26S)), get_total_genes_in_subg("Bhyb26S")), row.names=c("lost", "total"))
chisq.test(contingency_table)

#Just curious
contingency_table <- data.frame("D"=c(sum(unlist(abs_Bhyb26D)),sum(unlist(conserved))), "S"=c(sum(unlist(abs_Bhyb26S)), sum(unlist(conserved))), row.names=c("lost", "total"))
chisq.test(contingency_table)


#ABR113
contingency_table <- data.frame("D"=c(sum(unlist(abs_ABR113D)), get_total_genes_in_subg("ABR113D")), "S"=c(sum(unlist(abs_ABR113D)), get_total_genes_in_subg("ABR113S")), row.names=c("lost", "total"))
chisq.test(contingency_table)

#diploids
contingency_table <- data.frame("D"=c(sum(unlist(abs_Bd21)), get_total_genes_in_subg("Bdistachyon")), "S"=c(sum(unlist(abs_Bstacei)), get_total_genes_in_subg("Bstacei")), row.names=c("lost", "total"))
chisq.test(contingency_table)
##################################################################




############## Create lostgenes file for pseudogene pipeline ############## 

Ddf <- as.data.frame(unlist(abs_Bhyb26D))
Ddf <- rownames(Ddf)[which(Ddf[,1]==TRUE)]
Sdf <- as.data.frame(unlist(abs_Bhyb26S))
Sdf <- rownames(Sdf)[which(Sdf[,1]==TRUE)]
Ddf <- data.frame(synOg=Ddf, missing_subgenome="Bhyb26D", stringsAsFactors = F)
Sdf <- data.frame(synOg=Sdf, missing_subgenome="Bhyb26S", stringsAsFactors = F)
final <- rbind(Ddf, Sdf)
write.table(final, file="lostgenes.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
##################################################################


############## Create conservedgenes file for pseudogene pipeline ############## 

#Here we can see how many orthogroups fall into each category
sum(unlist(conserved))
sum(unlist(abs_Bhyb26D))
sum(unlist(abs_Bhyb26S))
sum(unlist(abs_ABR113D))
sum(unlist(abs_ABR113S))
sum(unlist(abs_Bd21))
sum(unlist(abs_Bstacei))
sum(unlist(abs_rice))


alldf <- as.data.frame(unlist(conserved))
allvec <- rownames(alldf)[which(alldf[,1]==TRUE)]
DTslim <- DT[which(DT$synOg %in% allvec)]

one_gene_per_genome <- function(mysynOg) {
  mini <- DTslim[synOg==mysynOg,]
  t <- table(mini$genome)
  if(sum(sapply(c('Bhyb26D', 'Bhyb26S', 'Bdistachyon', 'Bstacei'), function(g) unlist(unname(t[names(t)==g]))))==4) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#Takes about 6 minutes
tokeep <- sapply(allvec, function(s) one_gene_per_genome(s))
tokeep2 <- tokeep[which(tokeep==TRUE)]

get_minidfs <- function(mysynOg) {
  return( DTslim[synOg==mysynOg,] )
}
#Takes another 6 or 7 minutes
df_list <- lapply(names(tokeep2), function(s) get_minidfs(s))

gene_data <- lapply(df_list, function(m) {
  return(c(
    unique(m$synOg),
    unlist(unname(m[m$genome=="Bdistachyon", 'id'])),
    unlist(unname(m[m$genome=="Bhyb26D", 'id'])),
    unlist(unname(m[m$genome=="Bhyb26S", 'id'])),
    unlist(unname(m[m$genome=="Bstacei", 'id']))
  ))
})

final_gene_data <- do.call(rbind.data.frame, gene_data)
names(final_gene_data) <- c("synOg", "Bdgene", "Bhyb26Dgene", "Bhyb26Sgene", "Bsgene")
write.table(final_gene_data, file="conservedgenes.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")

save.image(file="/home/virginia/Desktop/bin_orthogroups.RData")

##################################################################



  



# 
# toPlot <- data.frame( 
#   Absent = factor(c('Bhyb26D', 'Bhyb26S', 'ABR113D', 'ABR113S', 'Bd21', 'ABR114'), levels=c('Bhyb26D', 'Bhyb26S', 'ABR113D', 'ABR113S', 'Bd21', 'ABR114')),
# #  Absent = factor(c('Bhyb26D', 'Bhyb26S', 'ABR113D', 'ABR113S', 'Bd21', 'ABR114', 'Osativa'), levels=c('Bhyb26D', 'Bhyb26S', 'ABR113D', 'ABR113S', 'Bd21', 'ABR114', 'Osativa')),
#   Value = c(
#             sum(unlist(abs_Bhyb26D)),
#             sum(unlist(abs_Bhyb26S)),
#             sum(unlist(abs_ABR113D)),
#             sum(unlist(abs_ABR113S)),
#             sum(unlist(abs_Bd21)),
#             sum(unlist(abs_Bstacei)) #,
#             #sum(unlist(abs_rice)))
# ))
# 
# 
# 
# 
# toPlot$mylabels <- ifelse(toPlot$Absent=='Bhyb26D', 'Bd-Bs-ABR113D-ABR113S------Bhyb26S-Os',
#         ifelse(toPlot$Absent=='Bhyb26S', 'Bd-Bs-ABR113D-ABR113S-Bhyb26D------Os',
#         ifelse(toPlot$Absent=='ABR113D', 'Bd-Bs------ABR113S-Bhyb26D-Bhyb26S-Os',
#         ifelse(toPlot$Absent=='ABR113S', 'Bd-Bs-ABR113D------Bhyb26D-Bhyb26S-Os',
#         ifelse(toPlot$Absent=='Bd21', '---Bs-ABR113D-ABR113S-Bhyb26D-Bhyb26S-Os',
#         ifelse(toPlot$Absent=='ABR114', 'Bd----ABR113D-ABR113S-Bhyb26D-Bhyb26S-Os',
#         NA  ))))))
# 
# toPlot$genome <-ifelse(toPlot$Absent=='Bhyb26D', 'Bhyb26',
#                 ifelse(toPlot$Absent=='Bhyb26S', 'Bhyb26',
#                 ifelse(toPlot$Absent=='ABR113D', 'ABR113',
#                 ifelse(toPlot$Absent=='ABR113S', 'ABR113',
#                 ifelse(toPlot$Absent=='Bd21', 'Bd21',
#                 ifelse(toPlot$Absent=='ABR114', 'ABR114',
#                 NA  ))))))
# 
# p <- ggplot(toPlot, aes(x = mylabels, y = Value, fill=genome)) +
#   geom_bar(stat = 'identity') +#fill="violetred")+#position = 'dodge') +
#   coord_flip()+
#   theme(
#     #axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.title.x = element_text(margin = margin(t = 10), size=15),
#     axis.text.x = element_text(size=14, angle=-70),
#     axis.text.y = element_text(size=13),
#     legend.text=element_text(size=14),
#     legend.title=element_text(size=15),
#     panel.background = element_blank(),
#     panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray90"),
#     panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
#     plot.margin=unit(c(1,1,1,1),"cm")) +
#   ylab("Number of Orthogroups") +
#   scale_y_continuous(expand = c(0, 0), breaks=seq(0, 450, by=50))+
#   scale_fill_manual("Missing genome", values=c('Bhyb26'='palevioletred1', 'ABR113'="slateblue1", 'Bd21'="violetred",  'ABR114'="goldenrod2")) #"dodgerblue4"))
# 
# ggsave(filename = "oct15_orthogroups.png", plot = p, width = 10, height = 4, dpi = 600, units = "in")
# 



















