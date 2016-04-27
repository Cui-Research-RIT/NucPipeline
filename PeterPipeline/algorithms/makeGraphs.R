

#Takes a LoDF, and makes a symmetrized graph for occupancy data
ocpGraph <- function(toGraph, xaxis, dataName, dataType, PNG_WIDTH=700, PNG_HEIGHT=1000){
	# Generate the filename
	png(paste0(dataName, "-", floor(length(xaxis)/2), ".png"), family="sans", width = PNG_WIDTH, height = PNG_HEIGHT)
	allColors <- c("red", "black", "blue", "green", "magenta", "orange", "darkgreen", "yellow")
	
	#Symmetrize
	toGraph <- (toGraph + toGraph[rev(rownames(toGraph)),])/2
	
	#Graph margins
	par(family="sans", mar=c(5, 6, 4, 2))
	chartmax <- max(toGraph)
	chartmin <- min(toGraph)
	# Start the plot with everything blank, add individual elements later
	if(dataType == "occupancy"){
		plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0.5, 3), col=rgb(0,0,0,0))
		mtext("Normalized Nucleosome occupancy", side=2, cex=3, line=3.2, family="sans")
		mtext(paste0(dataName, " Nucleosome profiles"), side=3, cex=3, line=1, family="sans")
	}else if(dataType == "dnase"){
		allColors <- c("red", "blue", "black")
		plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0, 18), col=rgb(0,0,0,0))
		mtext("Number of normalized sequence tags", side=2, cex=3, line=3.2, family="sans")
		mtext(paste0(dataName, " DNase-seq profiles"), side=3, cex=3, line=1, family="sans")
	}else if(dataType == "medip"){
		allColors <- c("red", "blue", "black")
		plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0, 5), col=rgb(0,0,0,0))
		mtext("Number of normalized sequence tags", side=2, cex=3, line=3.2, family="sans")
		mtext(paste0(dataName, " MeDIP profiles"), side=3, cex=3, line=1, family="sans")
	}
	
	# Axes
	axis(1, cex.axis=2.5)
	axis(2, cex.axis=2.5)
	# Labels
	mtext("Distance from p53 binding site (bp)", side=1, cex=3, line=3.6, family="sans")
	
	#Print all lines, one line per column
	for(colnum in 1:ncol(toGraph)){
		lines(xaxis, toGraph[,colnum], type="l", col=allColors[colnum], lwd=3)
		saveVec(toGraph[,colnum], paste0(dataType, "-", colnames(toGraph)[colnum]), dataName)
	}
	abline(v=0, col="black")
	# Legend of which datasets are which lines
	legend("topright", colnames(toGraph), lty=c(1,1), lwd=c(2.5,2.5),
		col=allColors[1:ncol(toGraph)], cex=2.5)
	dev.off()
}

# Takes a list with two data frames, which are the first 1000 rows of the 1000bp and 2000bp autocorrelation results
# and graphs the two data frames in to side-by-side graphs
corrGraph <- function(allData, dataName, PNG_WIDTH=700, PNG_HEIGHT=1000){
	
	allColors <- c("red", "black", "blue", "green", "magenta", "orange", "darkgreen", "yellow")
	xaxis <- 0:1000
	
	# Generate the filename
	
	
	# Keep track of iterations for scaling of various things based on 1000 or 2000 dataset
	iter <- 1
	
	# Plot the two graphs one at a time
	for(graphSet in allData){
		png(paste0(dataName, "-corr.png"), family="sans", width = PNG_WIDTH, height = PNG_HEIGHT)
	
		# Margins, two graphs
		par(family="sans", mar=c(5, 6, 5, 1)) 
		# Dynamically generate some lengths. This shouldn't change, but I'll leave it in in case we
		# decide to change the graphing window later
		len <- nrow(graphSet) - 1
		
		# Start the plot with everything blank, add individual elements later
		plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n",
			ylim=c(0,2500*iter), col=rgb(0,0,0,0))
		axis(1, cex.axis=2.5)
		axis(2, cex.axis=2.5)
		
#~ 		# Add the text specifying the window size
#~ 		if(iter == 1){
#~ 			mtext("40bp-avg", side=3, cex=2, adj=0, line=1, col="blue")
#~ 		}
		
		mtext("Distance from p53 binding site (bp)", side=1, cex=3, line=3.6, family="sans")
		mtext("Distance auto-correlation", side=2, cex=3, line=3.2, family="sans")
		mtext(paste0(dataName, "-", len*iter), side=3, cex=3, line=1, family="sans")
		abline(v=0, col="black")
		for(colnum in 1:(ncol(graphSet))){
			lines(xaxis, graphSet[,colnum], type="l", col=allColors[colnum], lwd=3)
			saveVec(graphSet[,colnum], paste0("distCorr", iter, "-", colnames(graphSet)[colnum]), dataName)
		}
		legend("topright", colnames(graphSet), lty=c(1,1), lwd=c(2.5,2.5),
			col=c(allColors[1:ncol(graphSet)], "black"), cex=2.5)
		
		for(nuc in 1:10){
			abline(v=(nuc * 100), col="grey")
		}
#~ 		abline(v=(which.max(p3K)+499), col="blue", lty=2)
		iter <- iter + 1
	}
	dev.off()
}

# This will generate the graphs for histone modification data from a LoDF
modGraph <- function(allData, xaxis, dataName, PNG_WIDTH=700, PNG_HEIGHT=1000){
	allColors <- c("red", "blue", "black", "green", "magenta", "orange", "darkgreen", "yellow")
	
	for(datacol in c(1, 2, 3)){
		datacols <- c(datacol, datacol + 3)
		graphSet <- allData[,datacols]
		legNames <- colnames(graphSet)
		
		datasetName <- strsplit(legNames[1], "-")[[1]][2]
		
		png(paste0(dataName, "-", datasetName, "-", xaxis[length(xaxis)], ".png"), family="sans", 
			width = PNG_WIDTH, height = PNG_HEIGHT)
		par(family="sans", mar=c(5, 6, 4, 1)) 
		#symmetrize
		graphSet <- (graphSet + graphSet[rev(rownames(graphSet)),])/2
		# Center the data on the x axis
		graphSet <- graphSet[xaxis + ceiling(nrow(graphSet)/2),]
		
		#Get the max of the datasets which are 
		chartmax <- max(graphSet)
		chartmin <- min(graphSet)
		if(datacol == 3){
			plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n",
				ylim=c(0, 50), col=rgb(0,0,0,0))
		}else if(datacol == 1){
			plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n",
				ylim=c(0, 1.5), col=rgb(0,0,0,0))
		}else if(datacol == 2){
			plot(xaxis, rep(0, length(xaxis)), type="l", xlab="", ylab="", xaxt="n", yaxt="n",
				ylim=c(1, 1.8), col=rgb(0,0,0,0))
		}
		axis(1, cex.axis=2.5)
		axis(2, cex.axis=2.5)
		mtext("Distance from p53 binding site (bp)", side=1, cex=3, line=3.6, family="sans")
		mtext("Number of normalized sequence tags", side=2, cex=3, line=3.2, family="sans")
		mtext(paste0(dataName, "-", xaxis[length(xaxis)]), side=3, cex=3, line=1, family="sans")
		
		for(colnum in 1:ncol(graphSet)){
			lines(xaxis, graphSet[,colnum], type="l", col=allColors[colnum], lwd=3)
			saveVec(graphSet[,colnum], paste0("modification-", colnames(graphSet)[colnum]), dataName)
		}
		
		legend("topright", legNames, lty=c(1,1), lwd=c(2.5,2.5),
			col=c(allColors[1:ncol(graphSet)], "black"), cex=2.5)
		
		abline(v=0, col="black")
		dev.off()
	}
}


tssBar <- function(tssData, dataName, PNG_WIDTH=500, PNG_HEIGHT=350){
	png(paste0(dataName, "-", "tssDistHistogram.png"), family="sans", width = PNG_WIDTH, height = PNG_HEIGHT)
	par(family="sans")
	barplot(tssData, space = 0, xaxt = 'n', yaxt='n', ylim = c(0, 0.5))
	mtext("Frequency", side=2, cex=2, line=2.2)
	mtext("Distance to TSS (kb)", side=1, cex=2, line=3)
	axis(side = 1, at = c(0:12), labels = c(NA, -25,-10,-5,-2,-1,0,1,2,5,10,25, NA), cex.axis=1.3)
	axis(side = 2, cex.axis=1.5)
	
	dataname18 <- paste0("../../tmp/", dataName, "18.txt")
	data18 <- read.delim(dataname18, header=F, as.is=T)
	len <- nrow(data18)
	mtext(paste0(dataName, " (", len, ")"), side=3, cex=2, line=1, family="sans")
	dev.off()
	
}



tssPie <- function(tssData, dataName, PNG_WIDTH=450, PNG_HEIGHT=450){
	png(paste0(dataName, "-", "tssDistPiechart.png"), family="sans", width = PNG_WIDTH, height = PNG_HEIGHT)
	par(family="sans")
	frac <- (tssData + rev(tssData))[1:6]
	names(frac) = c(">25 kb", "10-25 kb", "5-10 kb", "2-5 kb", "1-2 kb", "<1 kb")
	pie (frac, col = rainbow(6), radius = 1, cex = 1.2)
	
	dataname18 <- paste0("../../tmp/", dataName, "18.txt")
	data18 <- read.delim(dataname18, header=F, as.is=T)
	len <- nrow(data18)
	mtext(paste0(dataName, " (", len, ")"), side=3, cex=2, line=1, family="sans")
	dev.off()
	
}


saveVec <- function(dataRow, dataType, rowName){
	fname <- paste0("../", dataType, "-", length(dataRow))
	if(file.exists(fname)){
		prevData <- read.table(fname, sep=",", as.is=T, header=F)
		if(rowName %in% prevData[,1]){
			old <- which(prevData[,1] == rowName)
			prevData <- prevData[-old,]
		}
		dataRow <- c(rowName, dataRow)
		prevData <- rbind(prevData, dataRow)
		write.table(prevData, file=fname, quote=F, row.names=F, col.names=F, sep=",")
	}else{
		prevData <- data.frame(matrix(nrow=0, ncol=length(dataRow)))
		prevData <- rbind(prevData, dataRow)
		rownames(prevData)[1] <- rowName
		write.table(prevData, file=fname, quote=F, row.names=T, col.names=F, sep=",")
	}
}











