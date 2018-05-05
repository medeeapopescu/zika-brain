#!/bin/bash

###################################
# Using cufflinks to track fpkm
# of input files. The output file
# is further manipulated to remove 
# header and sort by first column
###################################

input_file=$1
output_file=$2
threads=$3
annotation=$4
scratch=$5

output_dir=cufflinks_tmp

echo $scratch
cd $scratch
date

# Loading necessary modules
module load use.singlecell
module load cufflinks/2.2.1.Linux_x86_64

# Using Cufflinks to compute fpkm
cufflinks \
--num-threads $threads \
--output-dir $output_dir \
--GTF $annotation \
--quiet $input_file

# delete the first row and sort by the first column (gene names)
sed 1d $output_dir/genes.fpkm_tracking | sort -k 1 > $output_file

echo
ls
echo 
ls $output_dir
echo

echo "Cufflinks Completed"
date
exit


