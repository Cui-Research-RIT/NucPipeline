source("../../algorithms/headTail.R")



bwigs <- dir()[strtail(dir(), 7) == ".bigWig"]

for(chrNum in c(1:22, "X")){
	chrNum <- paste0("chr", chrNum)
	for(bwig in bwigs){
		
		outname <- paste0("./outwigs/", chrNum, "-", strhead(bwig, -6), "wig")
		cmd <- paste0("./bigWigToWig ", bwig, " ", outname, " -chrom=", chrNum)
		print(cmd)
		system(cmd)
	}
}
