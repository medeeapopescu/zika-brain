#!/bin/bash

###################################
# Marking and removing duplicates
# using Piccard tools
###################################

input_file=$1
output_file=$2
output_read_prefix=$3
tool_dir=$4
metrics_file=$5
scratch=$6

echo $scratch
cd $scratch
echo
ls
echo
date
echo $input_file
echo $output_file
# Loading necessary modules
# module load use.singlecell
# module load samtools/1.1
source activate py3env

mkdir tmpdir

# mark duplicates with picard tools
# REMOVE_DUPLICATES=true do not write duplicates to the outputfile
# instead of writing duplicates to the output file with appropriate flags.
# default for strategy is SUM_OF_BASE_QUALITIES
# default for OPTICAL_DUPLICATE_PIXEL_DISTANCE is 100
# -Xmx16g is for allocating memory


picard  MarkDuplicates \
INPUT=$input_file \
OUTPUT=$output_file \
METRICS_FILE=$metrics_file \
\
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=true \
CREATE_INDEX=true \
TMP_DIR=tmpdir


# Using samtools to sort and index reads
# output of picard is already sorted
samtools index $output_file 

echo $output_read_prefix
# sorting by name -n and the last entry is output prefix
samtools sort -o $output_read_prefix -n $output_file 

echo
ls
echo
ls tmpdir
echo

echo "Marking Duplicates Completed"
date
exit


