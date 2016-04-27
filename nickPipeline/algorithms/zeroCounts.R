
# This function reads in a data file for a location with noncontinuous data, looks up its
# location in the location dataset, and fills in all gaps with 0s for -2000:2000 around the location
zeroProcess <- function(countfile, name18, name19){
	# Read in the name data
	hg18 <- read.delim(name18, header=F)
	hg19 <- read.delim(name19, header=F)
	
	filedata <- read.delim(countfile, header=F)
	posData <- hg18
	
	info <- strsplit(countfile, "=")[[1]]
	
	# Check if an hg19 dataset is being used
	if(strhead(info[2], 7) == "GM12878" || strhead(info[2], 4) == "K562"){
		posData <- hg19
	}
	
	#Get the location
	REpos <- posData[posData[,4] == info[1],3]
	
	fullData <- (REpos - 2000):(REpos + 2000)
	# Merge the vector generated above with the file data
	# This sets every missing value to NA
	fullData <- merge(fullData, filedata, by=1, all=T)
	# Set all NA values to 0
	fullData[is.na(fullData)] <- 0
	
	write.table(fullData, file=countfile, row.names=F, col.names=F, quote=F, sep="\t")
}


