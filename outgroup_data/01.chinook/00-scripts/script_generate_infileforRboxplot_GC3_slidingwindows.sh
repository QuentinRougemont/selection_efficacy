# TL - 130318
# bash script_generate_infileforRboxplot.sh Falbicollis.CDS.sumstats.sed.withgc3.withGC4fold.clean Falbicollis.CDS.sumstats.sed.withgc3.withGC4fold.clean.prboxplot.bidon 0 4000000
myfile=$1
outfile=$2
minimum=$(echo "$3" | bc)
slidingwindowslength=$(echo "$4" | bc)

#minimum=$(echo ".25") # first class, 0.25 sounds quite good

# initiate
rm tmp $outfile $outfile.genesID
start=$(echo "0")
end=$(echo "$minimum")

echo "GCclasses	nb_CDS	medianGC3	pNseq	pSseq	NbSynSites	totallengthCDS	pNpS" > $outfile

# adjust problems of the number of digits and then sort the file
grep -v "name" $myfile | awk '{print $11"	"$0}' | sed 's/^1	/1.0000 /g' | sed 's/^0	/0.0000  /g' | sed 's/^0.1	/0.1000	/g' | sed 's/^0.2	/^0.2000	/g' | sed 's/^0.3	/0.3000	/g' | sed 's/^0.4	/0.4000	/g' | sed 's/^0.5	/0.5000	/g' | sed 's/^0.6	/0.6000	/g' | sed 's/^0.7	/0.7000	/g' | sed 's/^0.8	/0.8000	/g' | sed 's/^0.9	/0.9000	/g' | sort -g -k 1,1 -t "	" > $myfile.sort

currentwindows=$(echo "1" | bc)
cutoffhalfslidingwindowslength=$(echo "$slidingwindowslength / 2" | bc)
totallengthprevious=$(echo "0" | bc) # initiate
medianGCseq=$(echo "NA")
switch=$(echo "0" | bc )
while read line; do 
	echo "$line" >> tmp
	currentlocus=$(echo "$line" | awk '{print $2}')
	lengthcurrent=$(echo "$line" | awk '{print $3}')
	gc3=$(echo "$line" | awk '{print $1}')
	totallengthcurrentwindows=$(echo "$totallengthprevious + $lengthcurrent" | bc )
	if [ "$totallengthcurrentwindows" -lt "$slidingwindowslength" ]; then
		echo "$currentwindows	$currentlocus	$gc3	$lengthcurrent	$totallengthcurrentwindows" >> $outfile.genesID
		totallengthprevious=$(echo "$totallengthcurrentwindows" | bc)
		if [ "$totallengthcurrentwindows" -ge "$cutoffhalfslidingwindowslength" ] && [ "$switch" -eq "0" ]; then # report the median GC value over the sliding windows
			medianGCseq=$(echo "$line" | awk '{print $12}')
			switch=$(echo "1" | bc )
		else
			continue
		fi
	else
		# compute stats
		nbgenes=$(less tmp | wc -l)
		pNseq=$(awk '{s+=$9} END {print s}' tmp | awk '{printf "%.4f\n", $1}') # $9 = pN
		pSseq=$(awk '{s+=$8} END {print s}' tmp | awk '{printf "%.4f\n", $1}') # $8 = pS
		NbSynSites=$(awk '{s+=$10} END {print s}' tmp | awk '{printf "%.1f\n", $1}') # $10 = Number of synonymous sites
		totallength=$(awk '{s+=$3} END {print s}' tmp | awk '{printf "%.1f\n", $1}') #$3 = total length of align
		pNpSratio=$(echo "scale=5; ($pNseq/($totallength-$NbSynSites))/($pSseq/$NbSynSites)" | bc)
		# print values
		echo "$currentwindows	$nbgenes	$medianGCseq	$pNseq	$pSseq	$NbSynSites	$totallength	$pNpSratio" >> $outfile
		# initiate again
		previouswindows=$(echo "$currentwindows")
		currentwindows=$(echo "$previouswindows + 1" | bc)
		totallengthcurrentwindows=$(echo "0" | bc) # reset
		totallengthprevious=$(echo "0" | bc) # reset
		switch=$(echo "0" | bc )
		medianGCseq=$(echo "NA")
		rm tmp
	fi
done < $myfile.sort
### print last
# generate the median GC value for the last one
totallength=$(awk '{s+=$3} END {print s}' tmp | awk '{printf "%.3f\n", $1}') #$3 = total length of align
cutoffhalftotallength=$(echo "$totallength / 2" | bc)
while read line; do
	if [ "$totallengthcurrentwindows" -ge "$cutoffhalftotallength" ] && [ "$switch" -eq "0" ]; then
		medianGCseq=$(echo "$line" | awk '{print $12}')
		switch=$(echo "1" | bc )
	else
		continue
	fi
done <  tmp
# compute stats
nbgenes=$(less tmp | wc -l)
pNseq=$(awk '{s+=$9} END {print s}' tmp | awk '{printf "%.3f\n", $1}') # $9 = pN
pSseq=$(awk '{s+=$8} END {print s}' tmp | awk '{printf "%.3f\n", $1}') # $8 = pS
NbSynSites=$(awk '{s+=$10} END {print s}' tmp | awk '{printf "%.3f\n", $1}') # $10 = Number of synonymous sites
pNpSratio=$(echo "scale=5; ($pNseq/($totallength-$NbSynSites))/($pSseq/$NbSynSites)" | bc)
# print values
echo "$currentwindows	$nbgenes	$medianGCseq	$pNseq	$pSseq	$NbSynSites	$totallength	$pNpSratio" >> $outfile
rm tmp $outfile.sort

# infile need to be like that:
#species Name    Size    N       S       P       W       Ps      Pn      NSS     D_Taj   GC3     scaffoldname    medianGC3sBRUT  medianGC3sNAexcluded    nbhaplotypeswithoutNA   name4fold       Size4fold       N4fold  S4fold  P4fold  W4fold  Ps4fold Pn4fold NSS4fold        D_Taj4fold      GC34fold
#Falbicollis_quantiles   NP_001268892.1.fst.clean.fst.clean.fst  612     40      2       0.449597        0.470196        0.0625  0.387097        154.25  -0.0827225      0.637255        NP_001268892.1  0.613   0.642   40      NP_001268892.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 285     40      0       0       0       0       0       100     -999    0.631579
#Falbicollis_quantiles   NP_001272809.1.fst.clean.fst.clean.fst  891     40      7       1.80011 1.64569 0.714754        1.08536 213.829 0.259826        0.484848        NP_001272809.1  0.444   0.463   40      NP_001272809.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 381     40      2       0.605372        0.470196        0.0974359       0.507937        129.75  0.542835        0.425197
#Falbicollis_quantiles   NP_001276784.1.fst.clean.fst.clean.fst  645     40      0       0       0       0       0       163.25  -999    0.548837        NP_001276784.1  0.5225  0.521   40      NP_001276784.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 279     40      0       0       0       0       0       98.75   -999    0.548387

