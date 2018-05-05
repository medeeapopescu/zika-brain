#!/bin/bash

## Rsync this example to the directory of your choice 
## To run, execute the following from the singlecell-submit head node:
## sbatch ./example_batch.sh

#SBATCH -N 1 ## modes
#SBATCH --ntasks=1             ## nodes  
#SBATCH --mem=256000 ## highest 
#SBATCH --cpus-per-task=12 ## cores, max of 12 on a node
#SBATCH -p bigmem                   ##bigmem: unresstricted;  general:2hr; long:24 hrs; unrestricted
#SBATCH --job-name=scde_frsh
#SBATCH --output="slurm-%A.out"      ## Output file (relative to current directory)  
 
srun R CMD BATCH scde_diff_exp_pairwise_epi.R
