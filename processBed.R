##A more optimal R script for
## processing a single bed file
## than originally used by peter for his pipeline

##is interactive only.

file<-readline(prompt="Filename ")
data<-read.table(file,header=TRUE)
drop<-c(4:10)
nd<-colnames(data)[drop]
#data[,4]<-NULL
data<-data[,!(names(data) %in% nd)]
data<-split(data,data[,1])
for(i in 1:length(data)){
        write.table(data[[i]], file=paste0(data[[i]][1,1],".txt"),sep="\t", row.names=F, col.names=F, quote=FALSE)
}

