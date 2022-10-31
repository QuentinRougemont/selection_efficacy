# TL - 130318
# bash script_generate_piNpiS.sh Berners_piNpiS.4fold.CDS.sumstats.final.clean
myfile=$1
#outfile=$2


#minimum=$(echo ".25") # first class, 0.25 sounds quite good

pNseq=$(awk '{s+=$8} END {print s}' $myfile | awk '{printf "%.1f\n", $1}')
pSseq=$(awk '{s+=$7} END {print s}' $myfile | awk '{printf "%.1f\n", $1}')
NbSynSites=$(awk '{s+=$9} END {print s}' $myfile | awk '{printf "%.1f\n", $1}')
totallength=$(awk '{s+=$2} END {print s}' $myfile | awk '{printf "%.1f\n", $1}')
pSratio=$(echo "scale=8; ($pSseq/$NbSynSites)" | bc)
pNratio=$(echo "scale=8; ($pNseq/($totallength-$NbSynSites))" | bc)
pNpSratio=$(echo "scale=6; ($pNseq/($totallength-$NbSynSites))/($pSseq/$NbSynSites)" | bc)

echo "$myfile	$pNseq	$pSseq	$NbSynSites	$totallength	$pSratio	$pNratio	$pNpSratio" 


# infile need to be like that:
#species Name    Size    N       S       P       W       Ps      Pn      NSS     D_Taj   GC3     scaffoldname    medianGC3sBRUT  medianGC3sNAexcluded    nbhaplotypeswithoutNA   name4fold       Size4fold       N4fold  S4fold  P4fold  W4fold  Ps4fold Pn4fold NSS4fold        D_Taj4fold      GC34fold
#Falbicollis_quantiles   NP_001268892.1.fst.clean.fst.clean.fst  612     40      2       0.449597        0.470196        0.0625  0.387097        154.25  -0.0827225      0.637255        NP_001268892.1  0.613   0.642   40      NP_001268892.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 285     40      0       0       0       0       0       100     -999    0.631579
#Falbicollis_quantiles   NP_001272809.1.fst.clean.fst.clean.fst  891     40      7       1.80011 1.64569 0.714754        1.08536 213.829 0.259826        0.484848        NP_001272809.1  0.444   0.463   40      NP_001272809.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 381     40      2       0.605372        0.470196        0.0974359       0.507937        129.75  0.542835        0.425197
#Falbicollis_quantiles   NP_001276784.1.fst.clean.fst.clean.fst  645     40      0       0       0       0       0       163.25  -999    0.548837        NP_001276784.1  0.5225  0.521   40      NP_001276784.1.fst.clean.fst.clean.fst.sites4foldonly.clean.fst 279     40      0       0       0       0       0       98.75   -999    0.548387

