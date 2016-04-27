#!/bin/bash -l

#SBATCH --mail-user jnf3769@rit.edu
#SBATCH --mail-type=FAIL
#SBATCH -t 185:0:0
#SBATCH -p work -n 1 -c 4
#SBATCH --mem=4000


Rscript pipeline.R occupancy $infolder 0 T F

