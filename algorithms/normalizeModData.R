
#This script was used once to normalize the modification datasets
#It should not need to be used again. Probably.
library(parallel)
source("../headTail.R")
orig <- getwd()
setwd("/mnt/TBSSD/nsmPipeline/datasets/modification/hg19")



datasets <- rep(0, 6)
names(datasets) <- c("GM12878-H3k4me3", "GM12878-H3k27me3", 
	"GM12878-H3k36me3", "K562-H3k4me3", "K562-H3k27me3", "K562-H3k36me3")

datasetsCt <- rep(0, 6)
names(datasetsCt) <- c("GM12878-H3k4me3", "GM12878-H3k27me3", 
	"GM12878-H3k36me3", "K562-H3k4me3", "K562-H3k27me3", "K562-H3k36me3")


inFiles <- dir()
inFiles <- inFiles[strtail(inFiles, 14) == "raw_signal.wig"]

for(inFile in inFiles){
	print(inFile)
	inData <- scan(inFile)
	nameData <- strsplit(inFile, "_")[[1]]
	dataset <- nameData[2]
	
	fSum <- sum(inData[c(F,T)])
	
	if(!dataset %in% names(datasets)){
		setwd(orig)
		stop("dataset doesn't match")
	}
	
	datasets[dataset] <- datasets[dataset] + fSum
	datasetsCt[dataset] <- datasetsCt[dataset] + 1
}

# Divide by hg19 genome length
# Got length from here: http://www.ncbi.nlm.nih.gov/assembly/2758/
#datasets <- datasets/3137144693

for(datNm in names(datasets)){
	datasets[datNm] <- datasets[datNm]/datasetsCt[datNm]
}

normFile <- function(inFile){
	print(inFile)
	inData <- scan(inFile)
	inData <- matrix(inData, ncol=2, byrow=T)
	nameData <- strsplit(inFile, "_")[[1]]
	dataset <- nameData[2]
	
	inData[,2] <- inData[,2]/datasets[dataset]*20
	
	outFile <- paste(nameData[1], nameData[2], "normalized.txt", sep="_")
	
	write.table(inData, file=outFile, row.names=F, col.names=F, quote=F, sep="\t")
}

mclapply(inFiles, normFile, mc.cores = getOption("mc.cores", 6L))
