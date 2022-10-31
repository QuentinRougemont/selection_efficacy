#!/bin/bash

#from the folder 37.PiPS_GC3s run me:
#exemple: ./01-scripts/01_run_piNpiS windowsize

size=$1 #window size

folder=03-piNpiS_$size
mkdir $folder

#original: for i in *sumstats.final ; 
for i in *sumstats ; 
do 
    name=$( basename $i )
    ./00-scripts/script_generate_infileforRboxplot_GC3_slidingwindows.sh $i \
    ${folder}/"${name}"."$size"_mb \
    0 \
    $size ; 
done 
