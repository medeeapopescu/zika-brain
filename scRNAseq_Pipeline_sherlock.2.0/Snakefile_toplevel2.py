# This snakemake file performs the following actions on RNAseq fastq data.
# Rarefaction curves for all the samples.
# All the samples are downsampled to the same depth.
# 
# This is how to call snakemake
#
# module load use.singlecell
# module load python/3.4.2
# snakemake -n -s snakefile argument (for running just the names)
# snakemake --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} \
#           --partition={params.partition} --mem={params.mem} -o slurm_output/%j-{params.name}" -p -j -k
# to kill snakemake process call killall snakemake
#
# Cool Commands
# nohup: doesn't kill a process if a shell command is closed
# git add -A: to update changed, added, and deleted file before committing
# use killall snakemake to kill processes you do not have access to
# 
# Revision History
# 2015.07.21 Brian Yu Created

# Import packages
import os, glob, subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
from snakemake.utils import read_job_properties

# Importing variables NEED TO CHANGE THIS ARGUMENT
root_folder = config["location"] 

# Importing relavent bio/sub-sample folders, IDs and variables
sample_table = pd.read_table(root_folder+"/code_analysis/cells.txt", header=0, index_col=0)
subsampleIDs = list(set(sample_table.index))
parameters = pd.read_table(root_folder+"/code_analysis/parameter.txt", index_col=0)
# print(parameters)

# Pulling out variables from parameter.txt
biosample = parameters.ix["biosample_name",'entry']
code_dir = parameters.ix["code_directory",'entry']
tool_dir = parameters.ix['tool_directory','entry']

# creating a list of number of reads to subsample
max_reads = int(parameters.ix['Down_Sample_Read_Number','entry'])
downsampleList = ["%5.0f" % ((x+1)/20*max_reads) for x in range(0,20)]
# print(downsampleList)

# Add include files or other snakefile rule files
include: "Snakefile.utils_Mark"
include: "Snakefile.utils_Felix"
include: "Snakefile_helper_Brian.py"
include: "Snakefile_import.py"
include: "Snakefile_singleCellAlignment.py"
include: "Snakefile_combineAlignmentResults.py"
include: "Snakefile_rarefactionCombine.py"

# User defined constants
workdir: "/local10G/wanxinw/rnaseq_results/"+parameters.ix["biosample_name",0]

#################################
# A list of all the rules
#################################

rule all:
  # sample should be a list of the cell barcode  names in a sequencing run. 
  # These are the long names in Miseq runs but ILxxxx-N7xx-N5xx in Nextseq runs
  input: 
    # Combined output
    "Combined_Analysis/rarefaction_gene_numbers.csv",
    "Combined_Analysis/rarefaction_read_proportion.2000000.csv"
    # Per cell output
    # expand("{subsample}/rseqc.{subsample}.geneBodyCoverage.txt", subsample=subsampleIDs),
    # expand("{subsample}/htseq.{subsample}.sortedByCoord.out", subsample=subsampleIDs),
    # expand("{subsample}/cufflinks.{subsample}.genes_fpkm_tracking.out", subsample=subsampleIDs)
    # fastqc of fastq reads in subsamples
    # expand("{subsample}/P1.{subsample}.fastqc_results.txt", subsample=subsampleIDs),
    # expand("{subsample}/P2.{subsample}.fastqc_results.txt", subsample=subsampleIDs)
  params: 
    name="RNAseq_top", 
    partition="quake,owners",
    qos="normal",
    time="5:00:00",
    mem="3000" 
  threads: 1
  version: "1.0"

