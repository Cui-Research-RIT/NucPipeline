

# This function will read in all location files in a folder and paste them together
# And then sort them by location. The second and third columns are used to compute 
# the centerpoint, and the third column is set to this centerpoint. The second column
# is set to the centerpoint - 1 for liftOver compatibility
generateLocs <- function(datafolder, assembly){
	orig <- getwd()
	setwd(datafolder)
	
	#Read in data
	locFiles <- dir()
	locFile <- data.frame(matrix(nrow=0, ncol=4))
	for(fl in locFiles){
		flData <- read.delim(fl, header=F)
		locFile <- rbind(locFile, flData)
	}
	
	setwd(orig)
	
	#Ensure all data is chr 1-22 and X
	chrs <- paste0("chr", 1:22)
	chrs <- c(chrs, "chrX")
	locFile[locFile[,1] %in% chrs,]
	
	#Use rowname as ID
	locFile[,4] <- rownames(locFile)
	
	
	dataName <- tail(strsplit(datafolder, "/")[[1]], n=1)
	#Save the new file
	outname <- paste0("./tmp/", dataName, assembly, ".txt")
	print(paste0("Saving data to ", outname))
	write.table(locFile, file=outname, quote=F, sep="\t", row.names=F, col.names=F)
	
}
