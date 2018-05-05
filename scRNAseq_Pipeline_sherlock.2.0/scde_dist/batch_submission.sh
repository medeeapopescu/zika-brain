#!/bin/bash

## Rsync this example to the directory of your choice 
## To run, execute the following from the singlecell-submit head node:
## sbatch ./example_batch.sh

#SBATCH -N 1 ## modes
#SBATCH --ntasks=1             ## nodes  
#SBATCH --mem=64000 ## highest 
#SBATCH --cpus-per-task=12 ## cores, max of 12 on a node
#SBATCH -p long                   ##bigmem: unresstricted;  general:2hr; long:24 hrs; unrestricted
#SBATCH --job-name=scde
#SBATCH --output="slurm-%A.out"      ## Output file (relative to current directory)  
 
srun R CMD BATCH --no-save --no-restore '--args input="/local10G/wanxinw/scde_job_submission/ERA_batch4_run1.2.3.4/for_tsne_batch_044.045.110.csv"' scde.submission.R batch_044.045.110.out
