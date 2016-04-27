
library(parallel)

setwd("hg19")

wigs <- dir()

procFl <- function(wig){
	
	x <- scan(wig, what="character", sep="\n")
	
	x <- as.numeric(x)
	x <- x[!is.na(x)]
	
	y <- 1:length(x)
	
	
	write.table(matrix(c(y,x),ncol=2), file=paste0(wig, "-raw"), row.names=F, col.names=F, quote=F, sep="\t")
}


mclapply(wigs, procFl, mc.cores = getOption("mc.cores", 3L))
