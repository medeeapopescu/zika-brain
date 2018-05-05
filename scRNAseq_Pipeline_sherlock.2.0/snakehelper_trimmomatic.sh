#!/bin/bash

###################################
# trim reads
###################################

read1=$1 # This will not be a .gz file
out.p1=$2
out.u1=$3
out.p2=$4
out.u2=$5
scratch=$6

source activate py3env
