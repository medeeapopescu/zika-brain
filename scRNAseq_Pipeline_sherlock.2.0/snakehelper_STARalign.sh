#!/bin/bash

###################################
# STAR alignment of 2 .gz files to
# a reference genome. Both names are
# passed in. Also passed in are other
# parameters as well as an output 
# directory
#### update for sherlock 2017.03.16
# cleaned up STAR and samtool executable directoreis by
# "source activate py3env"
###################################

P1=$1 # This will not be a .gz file
P2=$2 # This will not be a .gz file
genomedir=$3
threads=$4
outdir=$5
scratch=$6
#STAR=/local10G/resources/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR
source activate py3env

echo $scratch
cd $scratch
date

pwd

wc -l $P1
wc -l $P2

# STAR needs you to make this directory first
mkdir $outdir
echo $outdir

# STAR alignment
STAR \
--genomeDir $genomedir \
--readFilesIn $P1 $P2 \
--outSAMstrandField intronMotif \
--outFilterMultimapNmax 20\
--readFilesCommand cat \
--outSAMunmapped Within \
--outFileNamePrefix $outdir/ \
--runThreadN $threads \
--outReadsUnmapped Fastx \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate

# Index with samtools, default uses -b option and creates a .bam.bai file
samtools index $outdir/Aligned.sortedByCoord.out.bam

echo
ls
echo
ls $outdir
echo 

echo "STAR Alignment Completed"
date
exit


