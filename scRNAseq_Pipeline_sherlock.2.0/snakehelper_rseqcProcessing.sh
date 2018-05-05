#!/bin/bash

###################################
# RSeQC processing for quality check
# of RNAseq results
# 2017.04.10 Wanxin editted: RPKM_count.py -> FPKM_count.py as the former became obsolete per rseqc manual
###################################
input_file=$1
bed_file=$2
read_distribution_output=$3
FPKM_prefix=$4
RNA_fragment_size_output=$5
geneBody_prefix=$6
junction_prefix=$7
scratch=$8

echo $scratch
cd $scratch
date

# Currently must use wenying's version of RSeQC scripts
# source /local10G/wenyingp/resources/py2.7.9/bin/activate &&\
source activate py2env

# first find all required preix
echo "FPKM prefix is "$FPKM_prefix
echo "Gene Body prefix is "$geneBody_prefix
echo "Junction Annotation prefix is "$junction_prefix

read_distribution.py --input-file=$input_file --refgene=$bed_file > $read_distribution_output
echo "read distribution completed"
date

FPKM_count.py --input-file=$input_file --refgene=$bed_file --out-prefix=$FPKM_prefix > temp
echo
ls
echo
#mv $FPKM_prefix"_read_count.xls" $FPKM_prefix".xls"
echo "FPKM count completed"
date

RNA_fragment_size.py --input=$input_file --refgene=$bed_file > $RNA_fragment_size_output
echo "RNA fragment size completed"
date

geneBody_coverage.py --input=$input_file --refgene=$bed_file --format=png --out-prefix=$geneBody_prefix > temp
echo "geneBody completed"
date

junction_annotation.py --input-file=$input_file --refgene=$bed_file --out-prefix=$junction_prefix > temp
echo "junction annotation completed"
date

echo "RSeQC processing completed"
ls

exit


