###An R script for pulling down the 
##newer ENCODE data. Note that it requires a formal 
##list of filenames. 
###Can be invoked with Rscript
main<-function(args){
	fi<-args[1]
	fil<-file(fi)
	data2<-NULL
	data<-readLines(fil)
	for(item in data){
		#print(item)
		toGet<-paste0("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/",item)
		#print(toGet)
		#toGet<-gsub("[[:space:]]", "", toGet)
		
		#toDo<-paste0(c("wget ",toGet))
		#system(cat(toDo))
		
		data2<-c(data2,toGet)
	}
	close(fil)
	print(data2)
	write(data2, "encodeHTMLSforwgetHESC.txt")
}
args <- commandArgs(TRUE)
do.call(main, as.list(args))

