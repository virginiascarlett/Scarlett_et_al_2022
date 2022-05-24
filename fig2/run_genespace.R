library(GENESPACE)
gpar <- init_genespace(
  genomeIDs = c("Bdistachyon","Bstacei", "Bhyb26S", "Bhyb26D", "Osativa", "ABR113D", "ABR113S", "Bd28D", "Bd28S"),
  speciesIDs = c("Bdistachyon", "Bstacei", "Bhyb26S", "Bhyb26D", "Osativa", "ABR113D", "ABR113S", "Bd28D", "Bd28S"),
  versionIDs = c("v3.2", "v1.1", "v2.1", "v2.1", "v7.0", "v1.1", "v1.1", "v1.1", "v1.1"),
  ploidy = c(1,1,1,1,1,1,1,1,1), #ploidy of 1 means haploid genome
  wd = "/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new",
  nCores = 6,
  path2mcscanx = "/home/virginia/MCScanX",
  rawGenomeDir = "/home/virginia/phytozome"
)



#check_orthofinderinstall("orthofinder")
#Query is the genome with the most genes
#blkSize: how many genes are needed to call an SB
#synBuff: synteny buffer e.g. 100 genes around the syntenic block


#Parse annotations
parse_phytozome(gsParam = gpar, 
                genomeIDs = c("Bdistachyon", "Bstacei"),
                overwrite = T)

parse_annotations( gsParam=gpar,
                   genomeIDs = c("Bhyb26S", "Bhyb26D", "Bd28D", "Bd28S"),
                   headerEntryIndex = 3,
                   headerStripText = "locus=",
                   overwrite = T)

parse_annotations( gsParam=gpar,
                   genomeIDs = c("Osativa", "ABR113D", "ABR113S"),
                   headerEntryIndex = 4,
                   headerStripText = "locus=",
                   overwrite = T)


#Set the desired synteny parameters
#To see the full parameterization, print the synteny parameters as print(gpar$params$synteny)
gpar <- set_syntenyParams(gsParam = gpar) #could change params here, e.g. synBuff=50 or whatever

#Print the command to run orthofinder in an environment with orthofinder in the path
gpar <- run_orthofinder(gsParam = gpar, overwrite = T)

#FOR THE LOVE OF GOD SAVE HERE
#save.image(file = "/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/all_genomes_new.RData")
#load("/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/all_genomes_new.RData")

#gpar <- annotate_gff(gsParam = gpar, overwrite = T)

#Parse the orthofinder run into syntenic orthogroups
#For convenience, the syntenic block coordinates are returned, but not needed for downstream analyses
blks <- synteny(gsParam = gpar)
#write.table(blks, file="/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes/synBlocks.csv", row.names = FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(blks$params$synteny, file="/home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes_new/synBlocks.csv", row.names = FALSE, col.names=TRUE, quote=FALSE, sep=",")
#gpar <- pull_synOGs(gsParam = gpar)

#Build a pan-genome annotation
#The synteny constrained hits are used to predict their positions against a chosen reference genome
#This is meant to be the most useful output, everything is an intermediate file meant mainly for internal use
#pg <- pangenome(gsParam = gpar, refGenome = "Bdistachyon")
#pg <- pangenome(gsParam = gpar, refGenome = "Osativa")
#^Gives me the following error:
#Error in pangenome(gsParam = g par, refGenome = "Bdistachyon") : 
#  refGenome Bdistachyon not one of the genomeIDs 
#> gpar$genomeIDs
#NULL

#query_pangenome(pg, refChrom="Bd1", startOrder=10, endOrder=20)





#devtools::install_github("jtlovell/GENESPACE", upgrade = T)
