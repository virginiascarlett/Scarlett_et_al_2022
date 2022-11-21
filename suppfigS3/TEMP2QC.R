
read_my_file <- function(filename) {
  d <- read.table(filename, sep="\t", header = FALSE, stringsAsFactors = FALSE, fill=TRUE)  
  return(d)
}
get_medianQ <- function(myfile) {
  dat <- data.frame(read_my_file(myfile))
  dat <- dat[which(dat$V1==">>Per sequence quality scores")+1:nrow(dat),]
  qual <- dat[1:min(which(dat$V1==">>END_MODULE"))-1,]
  q <- as.numeric(qual$V2)
  return(median(q))
}
get_depth <- function(myfile) {
  dat <- data.frame(read_my_file(myfile))
  return(dat[5,]$V2)
}
setwd("/home/virginia/Documents/school/vogelLab/notebook/2021/TEpolymorphisms/QC")
files <- list.files(pattern = "^fastqc_data_sample.*\\.txt$")
medianQ <- sapply(files, function(f) get_medianQ(f))
final <- as.data.frame(medianQ)
final$depth <- as.numeric(unlist(unname(sapply(files, function(f) get_depth(f)))))

write.table(final, file="QCinfo.tsv", quote=FALSE, sep="\t", row.names=TRUE)

