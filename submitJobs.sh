#!/bin/bash
jobfile="./pipeline.sh"

for infolder in ./sites/Bo*
do
	jobname=$( echo $infolder | sed 's/.\/sites\///g' )
	#export outfile;
	export infolder;

	sbatch --qos=free -o ./slurmFiles/${jobname}.out \
	-e ./slurmFiles/${jobname}.err -J ${jobname} $jobfile
done
