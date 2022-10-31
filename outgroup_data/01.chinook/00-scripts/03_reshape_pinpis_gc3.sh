#!/bin/bash

windowssize=$1  #size of the windows 
#exemple windowssize=1000000 #for 1mb windows

OUTFOLDER="04-plot"
if [ ! -d "$OUTFOLDER" ]
then 
    echo "creating out-dir"
    mkdir "$OUTFOLDER"
fi

for i in 03-piNpiS_$windowssize/*_mb ; do grep -v "classes" $i |sed "s/^/$(basename $i)\t/g" >> 04-plot/$(basename $i).txt ; done

grep "classes" 03-piNpiS_$windowssize/*.CDS.sumstats."$windowssize"_mb |cut -d ":" -f 2- |uniq |sed 's/^/POP\t/g' > 04-plot/header.tmp
cat 04-plot/header.tmp 04-plot/*CDS.sumstats."$windowssize"_mb.txt > 04-plot/pnps_gc3.$windowssize.txt
rm 04-plot/*tmp

sed -i "s/_withoutquantiles.CDS.sumstats."$windowssize"_mb//g" 04-plot/pnps_gc3.$windowssize.txt
