#!/bin/bash

# testing some scripts

#SBATCH --job-name=BullTroutMPVA_JAGS
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --account=norock
#SBATCH --time=2-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcline@usgs.gov
#SBATCH -o %j.out

module load gnu8/8.3.0 R/4.1.1 jags/4.3.0

prun Rscript MPVA_ReddsOnly_Yeti.R