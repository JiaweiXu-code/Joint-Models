#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem 9000
#SBATCH --output=./cluster-out/02-%a.out
#SBATCH --error=./cluster-err/02-%a.err
#SBATCH --array=1-1040

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./programs/dataset.R ./cluster-logs/JM-$SLURM_ARRAY_TASK_ID.Rout

## run SAS command
sas -work /dev/shm -noterminal ./programs/JM.sas -log "./cluster-logs/nlmixed-$SLURM_ARRAY_TASK_ID.log" -print "./cluster-out/nlmixed-$SLURM_ARRAY_TASK_ID.lst" -sysparm "$SLURM_ARRAY_TASK_ID"
sas -work /dev/shm -noterminal ./programs/SJM.sas -log "./cluster-logs/nlmixed_nore-$SLURM_ARRAY_TASK_ID.log" -print "./cluster-out/nlmixed_nore-$SLURM_ARRAY_TASK_ID.lst" -sysparm "$SLURM_ARRAY_TASK_ID"
