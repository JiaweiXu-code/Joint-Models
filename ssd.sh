#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mem 2000
#SBATCH --output=./cluster-out/01-%a.out
#SBATCH --error=./cluster-err/01%a.err
#SBATCH --array=1-1000

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./ssd.R ./cluster-logs/ssd-$SLURM_ARRAY_TASK_ID.Rout
