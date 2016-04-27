
#Example to run:
#  srun --mincpus=16 --mem=75000 Rscript pipeline.R <DatasetFolderName> 3

#PETER COMMENT
# This function oversees the accumulation of all ChIP locations from chromosomes
# and the liftOver of the locations to hg19 if from hg18, or vice versa
# Parameters: datafolder, the name of the folder containing the location data
#trying to fix the x-11 issue

#NICK COMMENTS: 
# Runs Module 0...
# Converts assemblies via compiled liftover program
generateLocsD <- function(datafolder, assembly="18"){
	
	print("Compiling location data..")
	# This file contains the function to stitch the locations together and sort them
	source("algorithms/processInputLocs.R")
	generateLocs(datafolder, assembly)
	
	print(paste0("Lifting hg", assembly, " data..."))
	dataName <- tail(strsplit(datafolder, "/")[[1]], n=1)
	# Generate the names of the output files relative to the liftOver binary
	nameIN <- paste0("../tmp/", dataName, assembly, ".txt")
	nameOUT <- ""
	nameUM <- paste0(" ../tmp/", dataName, "UNMAPPED.txt")
	liftCmd <- ""
	if(assembly == "18"){
		nameOUT <- paste0("../tmp/", dataName, "19.txt")
		liftCmd <- paste0("./liftOver ", nameIN, " hg18ToHg19.over.chain ", nameOUT, nameUM)
	}else if(assembly == "19"){
		nameOUT <- paste0("../tmp/", dataName, "18.txt")
		liftCmd <- paste0("./liftOver ", nameIN, " hg19ToHg18.over.chain ", nameOUT, nameUM)
	}
	
	
	# Generates the liftOver command and executes it
	
	setwd("algorithms")
	
	for(liftFile in c(nameIN)){
		locFile <- read.delim(liftFile, header=F)
		#Get centerpoint
		locFile[,3] <- floor((locFile[,3] + locFile[,2])/2)
		locFile[,2] <- locFile[,3] - 1
		#Sort
		locFile <- locFile[order(locFile[,3]),]
		write.table(locFile, file=liftFile, quote=F, sep="\t", row.names=F, col.names=F)
	}
	
	system(liftCmd)
	setwd("..")
	print("Done.")
}
#NICK'S COMMENTS
# This section is by far the most time consuming
#Runs module 1. Changed from Peter's version in 2 places,
#Allowing only the occupancy to be run-- he suggested the first
# but never found the problem; I changed it permanently here 
extractDataD <- function(dataName, updateData){
	setwd("algorithms")
	#this is change 1 to make occupancy only
 	#datasetNames <- c("occupancy", "occurrence", "modification", "dnase", "medip")
	datasetNames <- c("occupancy")
	
	#Oddly enough, pipeline still depends on this declaration
	#need to detangle this, as it doesn't contribute to the 
	#occupancy analysis.
	hg19only <- c("modification", "dnase", "medip")
	
	outDir <- paste0("../tmp/", dataName, "/")
	
	#User may lie--check to see if updateData is possible
	if(!file.exists(outDir)){
		updateData <- F
	}
	#get rid of the old data files if we are replacing them
	#makes things easier in R than overriding
	if(!updateData){
		if(file.exists(outDir)){
			unlink(outDir, recursive=T)
		}
		dir.create(outDir, showWarnings = FALSE)
	}
	
	name18 <- paste0("../tmp/", dataName, "18.txt")
	name19 <- paste0("../tmp/", dataName, "19.txt")
	
	#Skip data that is already done
	if(updateData){
		extantData <- dir(outDir)
		datasetNames <- datasetNames[!datasetNames %in% extantData]
	}
	
	# In this case, only will be the occupancy analysis
	for(outname in datasetNames){
		
		if(updateData){
			dir.create(paste0(outDir, outname))
		}else{
			dir.create(paste0(outDir, outname), showWarnings = FALSE)
			dir.create(paste0(outDir, outname), showWarnings = FALSE)
		}
		
		#BENCHMARK: Takes ~45% of runtime to run the following code
		print(paste0("Extracting ", outname, " data..."))
		
		
		#This is the dependency of the hg19only variable
			#Test removal here

		#generate the system call to the perl script
		if(!outname %in% hg19only){
			extractCmd18 <- paste0("perl extractData.pl ", name18, " ../datasets/", 
				outname, "/hg18/ 2000 ", paste0(outDir, outname), "/")
			system(extractCmd18)
			extractCmd19 <- paste0("perl extractData.pl ", name19, " ../datasets/", 
				outname, "/hg19/ 2000 ", paste0(outDir, outname), "/")
			system(extractCmd19)
		}else{
			extractCmd19 <- paste0("perl extractData.pl ", name19, " ../datasets/", 
				outname, "/hg19/ 2000 ", paste0(outDir, outname), "/")
			system(extractCmd19)
		}
	}
	#####PERL SCRIPT SUMMARY (extractdata.pl
	#It takes four arguments: the tab-delimited file with the locations,
	# the folder with the raw data to be extracted, 
	#the window to extract around each position, and the output folder.
	#It reads in the tab-delimited file and processes it to a hash of arrays, 
	#where the keys of the hash are the chromosomes, 
	#and the arrays are lists of locations in each chromosome.
	# (the processCSV function)
	#The hash is then passed to the processFolder function, 
	#along with the other original parameters
	#It loops over all the files in parallel and writes out
	#accumulated data around particular positions once the file
	#reaches an appropriate size. 
	#	####IF LOOKING TO DO BACKEND DATABASE, HERE IS BEST#
	#
	#	#10-20-15: This portion of code crashed the research cluster
	#	#Edits made to .pl script to help prevent this 10-23-15
	###########
	
	setwd("..")
	
	cat("\n")
	zeroData(dataName, datasetNames)
	
	print("All done.")
}
#As far as I can tell, zeroData mostly just reformats the 
#data, inserting 0s in missing or implicit datapoints
#BENCHMARKIING--takes about 25% of runtime in series
zeroData <- function(dataName, datasetNames){
	
	
	library(parallel)
	source("algorithms/headTail.R")
	source("algorithms/zeroCounts.R")
	source("algorithms/fixWigOut.R")
	
	origDir <- getwd()
	dataDir <- paste0("./tmp/", dataName, "/")
	
	name18 <- paste0("../../", dataName, "18.txt")
	name19 <- paste0("../../", dataName, "19.txt")
	
	#another hg19 operation that can be removed
	hg19only <- c("modification", "dnase", "medip")
	
	for(outname in datasetNames){
		print(paste0("Zeroing gaps in dataset: ", outname))
		setwd(paste0(dataDir, outname))
		fnames <- dir()[!dir() %in% list.dirs(full.names=F)]
		
		if(outname == "modification"){
			mclapply(fnames, fixWigOut, name19, mc.cores = getOption("mc.cores", 16L))
		}else{
			#lapply(fnames, zeroProcess, name18, name19)
			mclapply(fnames, zeroProcess, name18, name19, mc.cores = getOption("mc.cores", 16L))
		}
		
		setwd(origDir)
	}
	
	
}
###Has no place in the occupancy pipeline, but 
# naiive removal causes failure. need to find
# what depends on this....
autoCorrD <- function(dataName, updateData){
	#uses mc.cores param from system to parallelize...
	#would benefit from higher core use in main pipeline, 
	# but since this version shouldn't run this, it is 
	# better to keep core usage low to reduce context switching
	library(parallel)
	
	origDir <- getwd()
	setwd(paste0("./tmp/", dataName, "/occurrence/"))
	
	if(!updateData){
		unlink("corr", recursive=T)
		unlink("sized", recursive=T)
		inFiles <- dir()
		dir.create("sized")
		
		print("Preparing occurrence data for autocorrelation")
		
		prepData <- function(fname){
			file2000 <- read.table(fname, header=F)
			file1000 <- file2000[1001:3001,]
			write.table(file1000, file=paste0("./sized/", fname, "-1000"), quote=F, row.names=F, col.names=F, sep="\t")
			write.table(file2000, file=paste0("./sized/", fname, "-2000"), quote=F, row.names=F, col.names=F, sep="\t")
		}
		
		mclapply(inFiles, prepData, mc.cores = getOption("mc.cores", 16L))
		
		dir.create("corr")
		setwd("sized")
		inFiles <- dir()
		
		print("Auto-correlating occurrence data")
		
		pb <- txtProgressBar(min = 0, max = length(inFiles), style = 3)
		
		runAutoCorr <- function(fname){
			setTxtProgressBar(pb, which(inFiles == fname))
			corrCmd <- paste0("perl ", origDir, "/algorithms/cross_correlation.pl ", fname)
			system(corrCmd)
		}
		
		mclapply(inFiles, runAutoCorr, mc.cores = getOption("mc.cores", 16L))
		print("Auto-correlation finished")
	}else{
		print("Auto-correlation skipped")
	}
	setwd(origDir)
	
}
#generic read-in function calling the other R functions
# used in more than one module
readData <- function(dataName, split25="both"){
	
	origDir <- getwd()
	source("algorithms/makeLoDF.R")
	source("algorithms/headTail.R")
	
	setwd(paste0("./tmp/", dataName))
	
	
	#....I have never tried with any individual calll here....
	splitNames <- c()
	if(split25 == "both"){
		splitNames <- read.delim(paste0("../", dataName, "18.txt"), header=F)[,4]
	}else if(split25 == "over25"){
		splitNames <- read.delim(paste0("../", dataName, "-over25-18.txt"), header=F)[,4]
	}else if(split25 == "under25"){
		splitNames <- read.delim(paste0("../", dataName, "-under25-18.txt"), header=F)[,4]
	}
	
	allData <- list()
	
	for(outname in dir()){
	
		backUp <- ".."
		if(outname == "occurrence"){
			setwd(paste0(outname, "/corr/"))
			backUp <- "../.."
		}else{
			setwd(outname)
		}
		
		print(paste0("Reading in ", outname, " data"))
		allData[[outname]] <- makeLoDF(splitNames)
		
		
		setwd(backUp)
	}
	
	
	setwd(origDir)
	
	return(allData)
}

graphD <- function(allData, dataName){
	
	source("algorithms/movAvg.R")
	source("algorithms/makeGraphs.R")
	
	print("Generating graphs and saving to output folder")
	
	origDir <- getwd()
	setwd("output")
#####SHOULD BE OK, NO NEED TO REMOVE, OVERRIDE WILL HAPPEN
#~ 	unlink(dataName, recursive=T)
#~ 	dir.create(dataName)
#~ 	setwd(dataName)
	
	#all data in this case will have only nucleosome occupancy
	for(dataType in names(allData)){
		
		dir.create(dataType, showWarnings=F)
		setwd(dataType)
		
		#Delete previous files
		unlink(list.files(pattern=dataName))
		
		symCharts <- c("occupancy", "modification", "dnase", "medip")
		#PETER:Occurrence graphs are significantly different and require different setup
		
		#TODO: REMOVE OCCURENCE
		if(dataType == "occurrence"){
			
			typeData <- allData[[dataType]]
			graphData <- data.frame(matrix(nrow=1001, ncol=0))
			
			for(datasetNum in 1:length(names(typeData))){
				dataset <- typeData[[datasetNum]][1:1001,]
				graphCol <- rowMeans(dataset)
				graphCol <- getAvg(graphCol, 40)
				
				graphData[,datasetNum] <- graphCol
				colnames(graphData)[datasetNum] <- names(typeData)[datasetNum]
			}
			graphData <- graphData[,order(colnames(graphData))]
			
			corrGraph(list(graphData[,strtail(colnames(graphData), 4)==1000], graphData[,strtail(colnames(graphData), 4)==2000]), dataName)
			
			
		}else if(dataType %in% symCharts){
			
			typeData <- allData[[dataType]]
			
			
			graphData <- data.frame(matrix(nrow=nrow(typeData[[1]]), ncol=0))
			
			for(datasetNum in 1:length(names(typeData))){
				dataset <- typeData[[datasetNum]]
				graphCol <- rowMeans(dataset)
				
				graphCol <- getAvg(graphCol, 60)
				if(dataType == "modification"){
					graphCol <- getAvg(graphCol, 60)
				}
				
				graphData[,datasetNum] <- graphCol
				colnames(graphData)[datasetNum] <- names(typeData)[datasetNum]
			}
			graphData <- graphData[,order(colnames(graphData))]
			
			if(dataType == "modification"){
				modGraph(graphData, -500:500, dataName)
				modGraph(graphData, -1000:1000, dataName)
			}else{
			#if(dataType == "occupancy" || dataType == "dnase"){
				ocpGraph(graphData[1501:2501,], -500:500, dataName, dataType)
				ocpGraph(graphData[1001:3001,], -1000:1000, dataName, dataType)
			}
			
		}
		
		setwd("..")
	}
	setwd(origDir)
	
}


tssDistD <- function(dataName, assembly="18"){
	
	origDir <- getwd()
	
	source("algorithms/tssDist.R")
	source("algorithms/headTail.R")
	source("algorithms/makeGraphs.R")
	
	tssData <- calcTssDist(dataName, assembly)
	setwd("output")
	dir.create("tssDist", showWarnings=F)
	setwd("tssDist")
	unlink(list.files(pattern=dataName))
	tssBar(tssData, dataName)
	tssPie(tssData, dataName)
	
	setwd(origDir)
	return(tssData)
}

####MAIN METHOD
#DEFAULTs: start from module 0, do not do progressive analysis, override old data, use hg18
# All can be overriden from commandline args
nsmPipe <- function(datafolder, module = 0, progressive = F, updateData=T, assembly="18", split25=F){
	###BENCHMARKING TOOL
	#startTime<-proc.time()
	
	##READ IN ARGS, convert to appropriate data type
	module <- as.numeric(module)
	
	dataName <- tail(strsplit(datafolder, "/")[[1]], n=1)
	
	progressive <- as.logical(progressive)
	updateData <- as.logical(updateData)
	split25 <- as.logical(split25)
	
	# PETER:Module 0 will use liftOver to generate the
	if(module == 0){
	
		generateLocsD(datafolder, assembly)
		
		print("Finished with module 0")
		if(progressive){
			module <- module + 1
		}
	}
	#PETER:
	# Module 1 involves determining the read locations of the data in hg18 and hg19,
	# extracting the appropriate subsequences, and otherwise getting the data ready for general use
	if(module == 1){
		if(updateData){
			extractDataD(dataName, updateData=T)
		}else{
			extractDataD(dataName, updateData=F)
		}
		print("Finished with module 1")
		if(progressive){
			module <- module + 1
		}
	}
	
	#PETER: Module 2 involves preparing the data for graphing, such as with the autocorrelation function
	if(module == 2){
		
	#	autoCorrD(dataName, updateData)
		print("Finished with module 2")
		if(progressive){
			module <- module + 1
		}
	}
	
	#PETER: Module 3 reads in the data and graphs it.
	if(module == 3){
		if(!split25){
			allData <- readData(dataName, "both")
			graphD(allData, dataName)
		}else{
			for(splitGrp in c("over25", "under25")){
				
				allData <- readData(dataName, splitGrp)
				graphD(allData, paste0(dataName, "-", splitGrp))
			}
		}
		
		
		
		print("Finished with module 3")
		if(progressive){
			module <- module + 1
		}
	}
	
	#PETER: Module 4 finds the closest TSS to each ChIP site and
	# Creates a chart based on distances
	if(module == 4){
		
		tssDistD(dataName, assembly)
		
		print("Finished with module 4")
		if(progressive){
			module <- module + 1
		}
	}
	#BENCMARKIING TOOL
	#print(proc.time-startTime)
}


# Command-line support!
#args <- commandArgs(TRUE)
#args2<-as.list(args)
#do.call(nsmPipe, as.list(args))
