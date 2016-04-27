#!/bin/bash -l

#SBATCH --mail-user afs2842@rit.edu
#SBATCH --mail-type=ALL
#SBATCH -t 72:0:0
#SBATCH -p work -n 1 -c 1
#SBATCH --mem=5000


Rscript pipeline.R $infolder 0 T F

