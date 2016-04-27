library(parallel)
source("../../algorithms/headTail.R")

setwd("outwigs")

wigs <- dir()
runWig <- function(wig){
	dataFile <- scan(wig, what="character", sep="\t")
	dataFile <- dataFile[strhead(dataFile, 12) != "variableStep"]
	dataFile <- as.numeric(dataFile)
	dataFile <- matrix(dataFile, ncol=2, byrow=T)
	write.table(dataFile, file=paste0("../doneWigs/", wig), col.names=F, row.names=F)
}

mclapply(wigs, runWig, mc.cores = getOption("mc.cores", 6L))
