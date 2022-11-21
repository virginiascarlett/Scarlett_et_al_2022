#I ran GENESPACE with just the two Bhyb26 subgenomes.
#Now I am parsing the output to get 1:1 homeologs.

library(ggplot2)
library(data.table)

get_og_comp <- function(myog) {
  return(DT[og==myog]$genome)
}

dat <- read.csv('/home/virginia/Documents/school/vogelLab/GENESPACE_Bhyb26/results/gffWithOgs.txt', sep='\t', stringsAsFactors = FALSE, header = TRUE)
DT <- as.data.table(dat)     
#takes about a minute
ogkey <- lapply(unique(dat$og), function(myog) get_og_comp(myog))
names(ogkey) <- unique(dat$og)
#temp <- unlist(lapply(ogkey, function(x) identical(as.vector(x), c("Bhyb26S", "Bhyb26D")))) #apparently there are no orthogroups with S, D, 
#so I will just ignore this vector
temp2 <- unlist(lapply(ogkey, function(x) identical(x, c("Bhyb26D", "Bhyb26S"))))
tokeep <- names(ogkey[which(temp2)])
DTslim <- DT[which(DT$og %in% tokeep),c("genome", "id", "og")]
library(tidyr)
final <- spread(DTslim, genome, id)
write.table(final[,c(2,3)], file="/home/virginia/Documents/school/vogelLab/GENESPACE_Bhyb26/Bhyb26_homeologs.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")











