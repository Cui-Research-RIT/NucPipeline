




# Command-line support!
args <- commandArgs(TRUE)
print("1")
print(args)
#args2<-as.list(args)
toRun<-args[1]
print("2")
print(toRun)
comboFlag<-FALSE
if(toRun=="all"){
	source("pipeTypes/pipelineAll.R")
	print("3")
} else if(toRun=="occupancy"){
	source("pipeTypes/pipelineOccupancy.R")
	print("4")
} else if(toRun=="occurrence"){
	source("pipeTypes/pipelineOccurrence.R")
	print("5")
} else if(toRun=="modification"){
	source("pipeTypes/pipelineModification.R")
	print("6")
} else if(toRun=="dnase"){
	source("pipeTypes/pipelineDNAse.R")
	print("7")
} else if(toRun=="medip"){
	source("pipeTypes/pipelineMedip.R")
	print("8")
} else if(toRun=="combo"){
	print("9")
	numAnalyses<-as.integer(args2[2])
	if(is.na(numAnalyses)||numAnalyses<=0){
			warning("analysis option 'combo' requires a subsequent, positive integer argument")
			stop("Invalid int following option combo")
	}
	listOfAnalyses<-args2[3:3+numAnalyses]
	args<-args[c(1:2,(3+numAnalyses+1)):length(args)]
	comboFlag<-TRUE
	source("pipeTypes/pipelineAll.R")
}else{
	print("10")
	warning("Invalid analysis option")
	stop("Analysis options are: occupancy, occurrence, modification, dnase, medip, or combo <number of analysises> <list>")
}
if(comboFlag==TRUE){
	datasetNames<-listOfAnalysis
}
print("11")
print(args)
args<-args[-1]
print(as.list(args))
do.call(nsmPipe, as.list(args))