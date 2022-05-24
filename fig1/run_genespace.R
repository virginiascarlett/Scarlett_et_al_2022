library(GENESPACE)
gpar <- init_genespace(
  genomeIDs = c("Bdistachyon","Bstacei", "Bhyb26S", "Bhyb26D", "Osativa", "ABR113D", "ABR113S"),
  speciesIDs = c("Bdistachyon", "Bstacei", "Bhyb26S", "Bhyb26D", "Osativa", "ABR113D", "ABR113S"),
  versionIDs = c("v3.2", "v1.1", "v2.1", "v2.1", "v7.0", "v1.1", "v1.1"),
  ploidy = c(1,1,1,1,1,1,1), #ploidy of 1 means haploid genome
  wd = "/home/virginia/Documents/school/vogelLab/GENESPACE_fig1",
  orthofinderMethod = "fast",
  nCores = 6,
  path2mcscanx = "/home/virginia/MCScanX",
  rawGenomeDir = "/home/virginia/phytozome"
)


#check_orthofinderinstall("orthofinder")
#Query is the genome with the most genes
#blkSize: how many genes are needed to call an SB
#synBuff: synteny buffer e.g. 100 genes around the syntenic block


#Parse annotations
#I'm doing each of these commands separately and in this order so
#that they come out in the order I want in the riparian plot.
#I'm sure there's a way to reorder the chromosomes in plot_riparian()
#but I don't know what it is.

parse_phytozome(gsParam = gpar, 
                genomeIDs = c("Bdistachyon", "Bstacei"),
                overwrite = T)

parse_annotations( gsParam=gpar,
                   genomeIDs = c("ABR113D", "ABR113S", "Osativa"),
                   headerEntryIndex = 4,
                   headerStripText = "locus=",
                   overwrite = T)

parse_annotations( gsParam=gpar,
                   genomeIDs = c("Bhyb26D", "Bhyb26S"),
                   headerEntryIndex = 3,
                   headerStripText = "locus=",
                   overwrite = T)





#Set the desired synteny parameters
#To see the full parameterization, print the synteny parameters as print(gpar$params$synteny)
gpar <- set_syntenyParams(gsParam = gpar) #could change params here, e.g. synBuff=50 or whatever

#Print the command to run orthofinder in an environment with orthofinder in the path
gpar <- run_orthofinder(gsParam = gpar, overwrite = T)

#Parse the orthofinder run into syntenic orthogroups
#For convenience, the syntenic block coordinates are returned, but not needed for downstream analyses
blks <- synteny(gsParam = gpar)

#FOR THE LOVE OF GOD SAVE HERE
#save.image(file = "/home/virginia/Documents/school/vogelLab/GENESPACE_fig1/fig1.RData")

blks$genomes$genomeIDs <- c("Bdistachyon", "ABR113D", "Bhyb26D", "Bhyb26S", "ABR113S", "Bstacei", "Osativa")
plot_riparian(gpar, plotRegions=F, blackBg = FALSE, genomeIDs=blks$genomes$genomeIDs)


plot_riparian(gpar, blackBg = FALSE, genomeIDs=blks$genomes$genomeIDs, refGenome="Bstacei")





