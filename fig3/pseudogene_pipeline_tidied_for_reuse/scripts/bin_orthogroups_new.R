#This script has five outputs: 
#First, an upset plot showing how many syntenic orthogroups are missing one particular genome.
#Second, a chi-squared test result describing whether two (sub)genomes are significantly
#different in terms of how many genes they've lost (i.e. subgenome dominance).
#Third, a file called lostgenes.txt which will serve as input for my script genespace2intervalsNEW.py.
#lostgenes.txt contains all the syntenic orthogroups of interest and looks like this: 
#synOg   missing_subgenome
#1       Bhyb26D
#2       Bhyb26D
#3       Bhyb26D
#9       Bhyb26D
#11      Bhyb26D
#13      Bhyb26D
#Fourth, a file called conserved_gene_intervals.bed which will allow me to run
#my pseudogene pipeline on conserved genes rather than missing genes, as a control.
#The format of this file is: (polyploid chrom, start, end, geneID of diploid syntenic ortholog)
#NOTE TO SELF: CHANGE THIS TO NAME OF PP GENE SO THAT GENEIDs WILL BE UNIQUE?
# BhD1    10   798 Bradi1g01820
# BhD1  3820 10475 Bradi1g01830
# BhD1 13585 20254 Bradi1g01840
# BhD1 22900 32043 Bradi1g01847
#A tetraploid will have two copies of every fully conserved gene, so each fully 
#conserved syntenic orthogroup will appear twice in this file.
#Also note that missing_gene_intervals.bed will had the diploid gene ID in the 
#last column, but conserved_gene_intervals.bed had the polyploid gene ID.
#Fifth, a file called conservedgenes.txt. This will be used by make_aliases.py
#to determine which polyploid genes are associated with which diploid genes.
#Looks like:
# Bhyb26D Bhyb26S Bdistachyon     Bstacei
# BrBhyv21000001m.g       BrBhyv21072774m.g       Bradi1g01820    Brast02G381900
# BrBhyv21000002m.g       BrBhyv21072773m.g       Bradi1g01830    Brast02G381700
# BrBhyv21000008m.g       BrBhyv21072768m.g       Bradi1g01840    Brast02G381600
# BrBhyv21000011m.g       BrBhyv21072766m.g       Bradi1g01847    Brast02G381400
#If there is only one progenitor, this file will only have three columns.

#User has to choose two polyploid subgenomes of interest. The whole pipeline is designed around
#the assumption that the user is working with a tetraploid with two well-defined subgenomes (separate files).
#If, for example, your two subgenomes of interest are ABR113D and ABR113S, then the pipeline will
#ultimately give you pseudogenes and dN/dS values for all the "missing genes" in ABR113.
#If you are interested in multiple polyploids, e.g. ABR113 and Bhyb26, you will have to re-run the 
#pipeline for each polyploid.

library(ggplot2)
library(data.table)
library(UpSetR)

get_og_comp <- function(myog) {
  return(DT[synOG==myog]$genome)
}
get_total_genes_in_subg <- function(mygenome) {
  return(
    length(unique(
      DT[which(DT$genome==mygenome),]$id
    ))
  )
}

########### EDIT HERE########### 
setwd("/home/virginia/Documents/school/vogelLab/notebook/2022/lost_genes")
dat <- read.csv('/home/virginia/Documents/school/vogelLab/GENESPACE_test/results/gffWithOgs.txt', sep='\t', stringsAsFactors = FALSE, header = TRUE)
my_genomes <- c("ABR113D", "ABR113S", "Bhyb26D", "Bhyb26S", "Bstacei", "Bdistachyon", "Osativa")
#GOI = Genome Of Interest. The two subgenomes of the polyploid you're interested in. Make sure these exactly match the corresponding entry in my_genomes.
GOI1 <- "Bhyb26D"
GOI2 <- "Bhyb26S"
progenitors <- c("Bdistachyon", "Bstacei") #Put these in the order so the first progenitor corresponds to GOI1 and the second to GOI2. 
#If there is only one progenitor, list it TWICE. 
############################################ 

#takes several minutes
#for every syntenic orthogroup, get the composition of the orthogroup in terms of genomes represented
#e.g.
#$`21902`
#[1] "ABR113D"     "ABR113S"     "Bd28D"       "Bd28S"       "Bdistachyon" "Bhyb26S"     "Bstacei"     "Osativa"  
DT <- as.data.table(dat) 
ogkey <- lapply(unique(dat$synOG), function(myog) get_og_comp(myog))
names(ogkey) <- unique(dat$synOG)


#For each possible combination of genomes where one genome is absent:
#Iterate over all syntenic orthogroups and check whether the genomes represented in that orthogroup
#correspond to the current combination of genomes.
#synog_results contains n elements where n is the length of my_genomes. Each element corresponds to a 
#particular combination of genomes, and has true/false values for every syntenic orthogroup.
#Takes about a minute.
synog_results <- lapply(seq(length(my_genomes)), function(n) 
  lapply(ogkey, function(synog_genomes) 
    all(my_genomes[-n] %in% synog_genomes) & !(my_genomes[n] %in% synog_genomes)
  )
)

#sum(unlist()) counts how many true values there are in the boolean list
num_missing_genes <- sapply(seq(length(synog_results)), function(n) sum(unlist(synog_results[n])) )
names(num_missing_genes) <- sapply(seq(length(synog_results)), function(n) paste(my_genomes[-n], collapse ="&") )


############## Create UpSet plot ############## 

tiff("nov2_upset.tiff", units="in", width=7, height=5.25, res=300)
upset(fromExpression(num_missing_genes),
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

##################################################################





############## Test for (pairwise) significance of biased gene loss (subgenome dominance) ############## 

gen1index <- which(my_genomes==GOI1)
gen2index <- which(my_genomes==GOI2)

contingency_table <- data.frame(
  "genome1"=c(num_missing_genes[gen1index], get_total_genes_in_subg(GOI1)), 
  "genome2"=c(num_missing_genes[gen2index], get_total_genes_in_subg(GOI2)),
  row.names=c("lost", "total"))

chisq.test(contingency_table)

##################################################################




############## Create lostgenes file for pseudogene pipeline ############## 

gen1index <- which(my_genomes==GOI1)
gen2index <- which(my_genomes==GOI2)

#TRUE values from synog_results gives us a vector of synOGs of interest
gen1df <- as.data.frame(unlist(synog_results[gen1index]))
gen1synogs <- rownames(gen1df)[which(gen1df[,1]==TRUE)]

gen2df <- as.data.frame(unlist(synog_results[gen2index]))
gen2synogs <- rownames(gen2df)[which(gen2df[,1]==TRUE)]

gen1final <- data.frame(synOG=gen1synogs, missing_subgenome=GOI1, stringsAsFactors = F)
gen2final <- data.frame(synOG=gen2synogs, missing_subgenome=GOI2, stringsAsFactors = F)

final <- rbind(gen1final, gen2final)
write.table(final, file="lostgenes.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
##################################################################


############## Create conservedgenes file for pseudogene pipeline ############## 

#Repeat what we did for synog_results, but now we are just getting those synOGs that
#have at least one gene from every genome, i.e. the gene is fully conserved. 
#NOTE: I decided to limit our conserved syntenic orthogroups to those that have exactly one
#gene per genome.
conserved <-  lapply(ogkey, function(synog_genomes) all(my_genomes %in% synog_genomes) & length(synog_genomes)==length(my_genomes) )

#Just get the TRUEs
alldf <- as.data.frame(unlist(conserved))
allvec <- rownames(alldf)[which(alldf[,1]==TRUE)]

#Make a copy of our original gffWithOgs dataframe that only contains fully conserved orthogroups.
#This makes for faster searching.
DTslim <- DT[which(DT$synOG %in% allvec)]

#Takes 5+ minutes!
#Create a list of dataframes, where each dataframe is the rows of DT with that particular synOG
df_list <- lapply(allvec, function(mysynOG)  DTslim[synOG==mysynOG,] )


conserved_GOI1 <- sapply(df_list, function(smalldf) {
  return(data.frame(
    chr=unname(unlist(smalldf[ which(smalldf$genome==GOI1), 'chr' ])),
    start=unname(unlist(smalldf[ which(smalldf$genome==GOI1), 'start' ])),
    end=unname(unlist(smalldf[ which(smalldf$genome==GOI1), 'end' ])),
    id=unname(unlist(smalldf[ which(smalldf$genome==GOI1), 'id' ])),
    stringsAsFactors = F
  ))})
conserved_GOI1 <- as.data.frame(t(conserved_GOI1))
  
conserved_GOI2 <- sapply(df_list, function(smalldf) {
  return(data.frame(
    chr=unname(unlist(smalldf[ which(smalldf$genome==GOI2), 'chr' ])),
    start=unname(unlist(smalldf[ which(smalldf$genome==GOI2), 'start' ])),
    end=unname(unlist(smalldf[ which(smalldf$genome==GOI2), 'end' ])),
    id=unname(unlist(smalldf[ which(smalldf$genome==GOI2), 'id' ])),
    stringsAsFactors = F
  ))})
conserved_GOI2 <- as.data.frame(t(conserved_GOI2))

#Takes a minute
gene_data <- lapply(df_list, function(small_df) {
  sapply(c(GOI1, GOI2, unique(progenitors)), function(GOI) {
    unlist(unname(small_df[small_df$genome==GOI, 'id']))
  })})
gene_data <- data.frame(do.call(rbind, gene_data))

finalbed <- rbind(conserved_GOI1, conserved_GOI2)
finalbed[] <- lapply(finalbed, function(x) type.convert(as.character(x)))
write.table(gene_data, file="conservedgenes.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(finalbed, file="conserved_gene_intervals.bed", row.names = FALSE, col.names=FALSE, quote=FALSE, sep="\t")

save.image(file="bin_orthogroups.RData")
