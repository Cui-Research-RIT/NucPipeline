
calcTssDist <- function(dataName, assembly){
	
	#t <- proc.time()
	#dataName <- "p53-154REs"
	#dataName <- "Chang_2014"
	
	if(assembly == "18"){
		dataname18 <- paste0("./tmp/", dataName, "18.txt")
		data18 <- read.delim(dataname18, header=F, as.is=T)
		hg18 <- read.delim("./datasets/tssPos18.txt", header=F, as.is=T)
	}else if(assembly == "19"){
		dataname18 <- paste0("./tmp/", dataName, "19.txt")
		data18 <- read.delim(dataname18, header=F, as.is=T)
		hg18 <- read.delim("./datasets/tssPos19.txt", header=F, as.is=T)
	}
	
	#Rather than subtracting differently for +/-
	#Easier to subtract one way and multiply by -1 if it should have been flipped
	hg18[,2][hg18[,2] == "+"] <- 1
	hg18[,2][hg18[,2] == "-"] <- -1
	hg18[,2] <- as.numeric(hg18[,2])
	
	#Initialize the buckets for the return values
	tssData <- rep(0, 12)
	tssData25 <- rep(0, 12)
	
	
	lastData <- list()
	lastStart <- rep(1, 25)
	bestLine <- ""
	names(lastStart) <- paste0("chr", c(1:22, "X", "Y", "M"))
	
	under25 <- c()
	
	#154 REs have their own data file
	if(dataName == "p53-154REs"){
		REdist <- scan("./datasets/p53-154REs-dist.txt")
		buckets25 <- cut(REdist, c(-Inf, -25000, -10000, -5000, -2000, -1000, 0,
			1000, 2000, 5000, 10000, 25000, Inf), right=F, labels=F)
		buckets <- cut(REdist, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, labels=F)
		for(bucket in buckets){
			tssData[bucket] <- tssData[bucket] + 1
		}
		for(bucket in buckets25){
			tssData25[bucket] <- tssData25[bucket] + 1
		}
		fracTssData <- tssData/length(REdist)
		fracTssData25 <- tssData25/length(REdist)
		
#~ 	}else if(dataName == "Su_2015"){
#~ 		REdist <- scan("./datasets/Su-dist.txt")
#~ 		buckets <- cut(REdist, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, labels=F)
#~ 		for(bucket in buckets){
#~ 			tssData[bucket] <- tssData[bucket] + 1
#~ 		}
#~ 		fracTssData <- tssData/length(REdist)
	}else{
		#Check each ChIP site
		for(dataRN in 1:nrow(data18)){
			
			dataRow <- data18[dataRN,]
#~ 			cat("\n\n\n")
#~ 			print(dataRow)
#~ 			print(lastData)
#~ 			print(lastStart)
			
			chr <- as.character(dataRow[1])
			pos <- as.numeric(dataRow[3])
			
			bestdist <- Inf
			
			if(chr %in% names(lastData)){
				bestdist <- (pos - lastData[[chr]][2]) * lastData[[chr]][1]
			}
			newLast <- lastData
			refCount <- lastStart[chr]
			check <- F
			if(refCount >= nrow(hg18)){
				bucket <- cut(bestdist, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, labels=F)
				bucket25 <- cut(bestdist, c(-Inf, -25000, -10000, -5000, -2000, -1000, 0,
					1000, 2000, 5000, 10000, 25000, Inf), right=F, labels=F)
				tssData[bucket] <- tssData[bucket] + 1
				tssData25[bucket25] <- tssData25[bucket25] + 1
				
				break
			}
			
			#Both lists are sorted- no need to start at the beginning for each
			for(refRN in refCount:nrow(hg18)){
				refRow <- hg18[refRN,]
				
				refChr <- as.character(refRow[1])
				strand <- as.numeric(refRow[2])
				tss <- as.numeric(refRow[3])
				
				#Add the new data to the list
				refCount <- refCount + 1
				#New data was correct chromosome, check score
				if(chr == refChr){
					check <- T
					newLast[[refChr]] <- c(strand, tss)
				}
				
				#Most recent data was correct chromosome, check if better or not
				if(check){
					bestLine <- refRow
					check <- F
					newdist <- (pos - newLast[[chr]][2]) * newLast[[chr]][1]
					
					#Better score! Write out new starting points and keep going
					if(abs(newdist) < abs(bestdist)){
						
						bestdist <- newdist
						lastData <- newLast
						lastStart[chr] <- refCount
#~ 						print(bestLine)
					#Worse score- last best is the minimum. Write distance to buckets, and move on
					#to next site without overwriting starting points
					}else if(abs(newdist) > abs(bestdist)){
						if(bestdist < 25000){
							under25 <- c(under25, dataRN)
						}
						bucket <- cut(bestdist, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, labels=F)
						bucket25 <- cut(bestdist, c(-Inf, -25000, -10000, -5000, -2000, -1000, 0,
							1000, 2000, 5000, 10000, 25000, Inf), right=F, labels=F)
						tssData[bucket] <- tssData[bucket] + 1
						tssData25[bucket25] <- tssData25[bucket25] + 1
						
						break
					}
				}
				if(refCount >= nrow(hg18)){
					bucket <- cut(bestdist, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, labels=F)
					bucket25 <- cut(bestdist, c(-Inf, -25000, -10000, -5000, -2000, -1000, 0,
						1000, 2000, 5000, 10000, 25000, Inf), right=F, labels=F)
					tssData[bucket] <- tssData[bucket] + 1
					tssData25[bucket25] <- tssData25[bucket25] + 1
					break
				}
			}
			
		}
		fracTssData <- tssData25/nrow(data18)
	}
	
	# Write out the data into the csv file
	
	outRow <- paste0(tssData, " (", sprintf("%1.2f%%", 100*fracTssData), ")")
	outRow25 <- paste0(tssData25, " (", sprintf("%1.2f%%", 100*fracTssData), ")")
	dataTable <- "./output/tssPositions5.csv"
	if(file.exists(dataTable)){
		prevData <- read.table(dataTable, sep=",", as.is=T)
		if(dataName %in% prevData[,1]){
			old <- which(prevData[,1] == dataName)
			prevData <- prevData[-old,]
		}
		outRow <- c(dataName, outRow)
		prevData <- rbind(prevData, outRow)
		write.table(prevData, file=dataTable, quote=T, row.names=F, col.names=F, sep=",")
	}else{
		prevData <- data.frame(matrix(nrow=0, ncol=12))
		prevData <- rbind(prevData, outRow)
		colnames(prevData) <- levels(cut(0, c(-Inf, seq(-5000, 5000, by=1000), Inf), right=F, dig.lab=5))
		rownames(prevData)[1] <- dataName
		write.csv(prevData, file=dataTable, quote=T)
	}
	
	dataTable <- "./output/tssPositions25.csv"
	if(file.exists(dataTable)){
		prevData <- read.table(dataTable, sep=",", as.is=T)
		if(dataName %in% prevData[,1]){
			old <- which(prevData[,1] == dataName)
			prevData <- prevData[-old,]
		}
		outRow25 <- c(dataName, outRow25)
		prevData <- rbind(prevData, outRow25)
		write.table(prevData, file=dataTable, quote=T, row.names=F, col.names=F, sep=",")
	}else{
		prevData <- data.frame(matrix(nrow=0, ncol=12))
		prevData <- rbind(prevData, outRow25)
		colnames(prevData) <- levels(cut(0, c(-Inf, -25000, -10000, -5000, -2000, -1000, 0,
		1000, 2000, 5000, 10000, 25000, Inf), right=F, dig.lab=5))
		rownames(prevData)[1] <- dataName
		write.csv(prevData, file=dataTable, quote=T)
	}
	
	if(dataName != "p53-154REs"){
		dataU25 <- data18[under25,]
		dataO25 <- data18[-under25,]
		write.table(dataU25, file=paste0("./tmp/", dataName, "-under25-18.txt"), row.names=F, col.names=F, quote=F, sep="\t")
		write.table(dataO25, file=paste0("./tmp/", dataName, "-over25-18.txt"), row.names=F, col.names=F, quote=F, sep="\t")
	}
	
	#print(proc.time() - t)
	return(fracTssData)
}


