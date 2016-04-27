# This function takes a vector and a window size
# And takes a moving average acros the vector with that window
# The final vector of the same length is returned
getAvg <- function(inVec, win = 30){
	
	retData <- rep(0, length(inVec))
	for(i in 1:length(inVec)){
		low <- i - floor(win/2)
		high <- i + (floor(win/2) - 1)
		if(low < 1){low <- 1}
		if(high > length(inVec)){high <- length(inVec)}
		retData[i] <- mean(inVec[low:high])
	}
	return(retData)
}
