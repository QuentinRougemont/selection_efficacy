#Diploid
cut -f 1  CDS.start.end.txt|grep "NC" |uniq |head -n 30 > list_NC
mkdir DIPLOID
awk 'NR==FNR { pat[$0]=1 } NR>FNR { for (p in pat) if ($0 ~ p) {print;next} }' list_NC CDS.start.end.txt > DIPLOID/NC_CDS.stats.end.txt
#faster than awk:
 grep -Ff list_NC CDS.start.end.txt > DIPLOID/NC
 cut -f 4 DIPLOID/NC |uniq > DIPLOID/uniq.cds.NC.txt

#keep only CDS in chr1-30 for each pop
for i in $(cat DIPLOID/uniq.cds.NC.txt) ; do
    for j in $(cat listpop) ; do
         grep "$i" "$j"_withoutquantiles.CDS.sumstats >> DIPLOID/"$j"_withoutquantiles.CDS.sumstats.DIPLOID;
    done ; 
done

#faire pinpis global
cd DIPLOID
for i in *DIPLOID ; do ../00-scripts/00_script_generate_piNpiS.sh $i >>pnps.diploidonly.txt;done


#script to perform bootstrap analysis

mkdir 06.boot
cd 06.boot

## bootstrap:
for i in $(cat listpop) ; do
    mkdir $i
    cd $i
    echo $(pwd)
    Rscript ../boot.R "$i"_withoutquantiles.CDS.sumstats.DIPLOID $i ;
    cd ../
done

for i in */*.random*txt ; do ../00-scripts/00_script_generate_piNpiS.sh $i >> boot.pnps.txt ; done
echo "POP myfile   pNseq  pSseq  NbSynSites     totallength    pSratio        pNratio        pNpSratio" >> header
cat header boot.pnps.txt >> tmp
mv tmp boot.pnps.txt
sed -i "s/\//\t/g" boot.pnps.txt
