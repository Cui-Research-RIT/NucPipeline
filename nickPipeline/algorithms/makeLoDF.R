
#Remove the largest 2% of data, which is often outliers
remOutliers <- function(allData){
	for(dataset in names(allData)){
		skewData <- allData[[dataset]]
		skewData <- skewData[,order(skewData[400,])]
		cap <- ceiling(ncol(skewData)*0.98)
		skewData <- skewData[,1:cap]
		allData[[dataset]] <- skewData
	}
	return(allData)
}

# This function will read in all files in the current working directory
# and accumulate them into a list of data frames (LoDF) where the list index is the dataset
# and each data frame contains all of the data for the files of that dataset, one file per column
# The LoDF is returned
makeLoDF <- function(splitNames){
	
	inFiles <- dir()
	datasets <- c()
	chipLocs <- c()
	
	#Get names of datasets from the filenames
	for(datapt in strsplit(inFiles, "=")){
		chipLocs <- c(chipLocs, datapt[1])
		datasets <- c(datasets, datapt[2])
	}
	
	#Get the unique names of the datasets. Sort them for consistency
	datasets <- sort(unique(datasets))
	
	# Initialize the list with all dataset names
	allData <- list()
	for(dataset in datasets){
		if(dataset == "widom"){
			next
		}
		allData[[dataset]] <- list()
	}
	
	
	#print(c(address(allData), refs(allData)))
	pb <- txtProgressBar(min = 0, max = length(inFiles), style = 3)
	
	for(i in 1:length(inFiles)){
		setTxtProgressBar(pb, i)
		
		inFile <- inFiles[i]
		
		fileData <- strsplit(inFile, "=")[[1]]
		locName <- fileData[1]
		dataset <- fileData[2]
		if(dataset == "widom"){next}
		if(!locName %in% splitNames){next}
		
		# By reading in the file with scan and adding it to a list, R is more efficient and runs ~3x faster
		allData[[dataset]][[locName]] <- scan(inFile, quiet=T)[c(F,T)]
		
	}
	
	print("All files read. Converting list to data frame")
	# Convert the list of lists to a list of data frames for ease of use later
	for(dataset in names(allData)){
		allData[[dataset]] <- as.data.frame(allData[[dataset]])
	}
	allData <- remOutliers(allData)
	# Return the LoDF
	return(allData)
}
