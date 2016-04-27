
#Modification dataset is for every 20 bp
#This distributes the value along that length
fixWigOut <- function(countfile, name19){
	
	# Read in the name data
	posData <- read.delim(name19, header=F)
	
	filedata <- read.delim(countfile, header=F)
	
	info <- strsplit(countfile, "=")[[1]]
	
	#Get the location
	REpos <- posData[posData[,4] == info[1],3]
	
	fullDataPos <- (REpos - 2000):(REpos + 2000)
	# Merge the vector generated above with the file data
	# This sets every missing value to NA
	fullData <- merge(fullDataPos, filedata, by=1, all=T)
	# Set all NA values to 0
	fullData[is.na(fullData)] <- 0
	
	# Fill in window of 20
	for(wigRow in 1:nrow(filedata)){
		pos <- filedata[wigRow, 1]
		val <- filedata[wigRow, 2]
		if(pos %in% fullDataPos){
			for(newpos in (pos-10):(pos+9)){
				fullData[fullData[,1] == newpos,2] <- val
			}
		}
	}
	
	write.table(fullData, file=countfile, row.names=F, col.names=F, quote=F, sep="\t")
}

