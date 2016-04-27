###processBedLot is a significantly more useful 
#version of processBed.R. It is invoked via R script
#and takes a folder of bed files. It converts them
#to the pipeline acceptable version and appends a tag to 
#let you know it is done. This tag comes in handy later for
#renaming the files and trimming off unnecessary info


main<-function(file,ins){
	data<-read.table(file,header=TRUE)
	drop<-c(4:10)
	nd<-colnames(data)[drop]
	#data[,4]<-NULL
	data<-data[,!(names(data) %in% nd)]
	data<-split(data,data[,1])
	setwd(ins)
	for(i in 1:length(data)){
        	write.table(data[[i]], file=paste0(data[[i]][1,1],".txt"),sep="\t", row.names=F, col.names=F, quote=FALSE)
	}
}
args<-commandArgs()
print(args)
setwd(args[6])
flist<-list.files()
starting<-getwd()
for(f in 1:length(flist)){
	system(paste0("mkdir ", flist[f], "_processed"))
	ins<-paste0(flist[f],"_processed")
	main(flist[f],ins)
	setwd(starting)
	#system(paste0("mkdir ", flist[f], "_processed"))
	#system(paste0("mv *.txt ", flist[f],"_processed"))	
	#tomove<-list.files(pattern=".txt")
	#print(flist[f])
	#for(file in tomove){
	#	print(tomove)
		#print(paste0(getwd(),"/",flist[f],"_processed/",flist[f]))
		#file.rename(file,paste0(getwd(),"/",flist[f],"_processed/",flist[f]))
		#system(paste0("mv ", file, " ", "/",flist[f],"_processed"))
	#}
}
