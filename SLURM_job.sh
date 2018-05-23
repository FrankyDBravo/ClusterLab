#!/bin/bash
# Job name:
#SBATCH --job-name=test1
#
# Account:
#SBATCH --account=fc_evodyn
#
# Partition:
#SBATCH --partition=savio
#
# Tasks per node
#SBATCH --ntasks-per-node=20
#
# Nodes
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run:
module load gcc openmpi # or module load intel openmpi
ht_helper.sh -t taskfile -r 20