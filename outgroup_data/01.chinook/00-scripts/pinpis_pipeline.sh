

#pipeline to estimate piNpiS:

# first obtain the gff:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/021/735/GCF_002021735.2_Okis_V2/GCF_002021735.2_Okis_V2_genomic.gff.gz

# extract the gene:
zcat GCF_002021735.2_Okis_V2_genomic.gff.gz |\
   grep -v "#" |awk '$3=="gene" {print $0}' > gene.gff
zcat GCF_002021735.2_Okis_V2_genomic.gff.gz |\
   grep -v "#" |awk '$3=="CDS" {print $0}' > CDS.gff 
sed 's/;/\t/g' CDS.gff  |cut -f 1,3-5,9,10,11,12,14 |less 
# as we can see multiple CDS are repeated for a same gene. Why?

sed 's/;/\t/g' CDS.gff  |cut -f 1,3-5,9,10,11,12,14 |\
    awk '{print $1"_"$3"_"$4"\t"$0 }' |awk '!seen[$1]++' |wc -l

#this is problematic and we remove them with awk:
awk '!seen[$1,$4,$5]++' CDS.gff > CDS.reshaped.gff

#we now produce one vcf for each pop.
#A full vcf was produced with gatk pipeline, and filtered based on quality score
for i in 01.scripts/00.extract_pop*.sh ; do msub $i ; done 

#then we will simplifiy each vcf to keep only the cds from the vcf:

zcat GCF_002021735.2_Okis_V2_genomic.gff.gz |\
   grep -v "#" |awk '$3=="CDS" {print $0}' |\
   sed 's/;/\t/g' |cut -f 1,3-5,9,10,11,12,14 |less

#### reshape the gff:
#use thibault script


#check the results:
awk '$2 % 3 != 0 {print $0}' check_length.txt | uniq |wc -l



################## RESHAPE GFF V1 ######################################
grep -v "#" GCF_002021735.2_Okis_V2_genomic.gff |cut -f 1,3,4,5,7,8,9 |sed 's/;/\t/g' |cut -f 1,2,7,8 >> colomns.to_reshape
grep -v "#" GCF_002021735.2_Okis_V2_genomic.gff |cut -f 2-8 > col2_8_to.keep
cut -f 1 colomns.to_reshape |sed -e 's/\.2$//g' -e 's/\.1$//g' >  colon1
cut -f 4 colomns.to_reshape > col9_2
cut -f 3 colomns.to_reshape |sed 's/-[0-9]*$//g' > col9_1
grep -v "#" GCF_002021735.2_Okis_V2_genomic.gff |cut -f 9 |cut -d ";" -f 3- > col9_3
paste col9_1 col9_2 col9_3 |sed 's/\t/;/g' > col9_reshaped
paste colon1 col2_8 col9_reshaped >> gff.reshaped
paste colon1 col2_8 col9_reshaped >> gff.reshaped


mkdir reshaped2
cd reshaped2
zcat ../gff.reshaped.gz |cut -f 1-8 > col1_8
zcat ../gff.reshaped.gz |cut -f 9  |sed 's/\.1;/;/g' |sed 's/\.2;/;/g' > newcol9
paste col1_8 newcol9 > gff.reshaped.2



zcat gff.reshaped2.gz |cut -f 1,3,4,5,9 |cut -d ";" -f 1 |grep "CDS" |awk '{if (seen[$1,$2,$3,$4]++) print $0,"==>","failed"; else print $0,"=>","Pass"; }' |less

zcat gff.reshaped2.gz |cut -f 1,3,4,5,9 |cut -d ";" -f 1 |grep "CDS" |awk '{if (seen[$1,$2,$3,$4]++) print $0,"==>","failed"; else print $0,"=>","Pass"; }' |grep "faile" |cut -f 5 |sort |uniq |wc -l
zcat gff.reshaped2.gz |cut -f 1,3,4,5,9 |cut -d ";" -f 1 |grep "CDS" |awk '{if (seen[$1,$2,$3,$4]++) print $0,"==>","failed"; else print $0,"=>","Pass"; }' |grep "Pass" |cut -f 5 |sort |uniq |wc -l


zcat gff.reshaped2.gz |cut -f 1,3,4,5,9 |cut -d ";" -f 1 |grep "CDS" |awk '{if (seen[$1,$2,$3,$4]++) print $0,"==>","failed"; else print $0,"=>","Pass"; }' |grep "Pass" |cut -f 5 |sort |uniq > CDS.pass.txt
zcat gff.reshaped2.gz |cut -f 1,3,4,5,9 |cut -d ";" -f 1 |grep "CDS" |awk '{if (seen[$1,$2,$3,$4]++) print $0,"==>","failed"; else print $0,"=>","Pass"; }' |grep "fail" |cut -f 5 |sort |uniq > CDS.failed.txt

CAP, number of CDS: 89265
Total number of CDS: 89365
Number of CDS non-duplicated: 


rsync -avh TSO/wanted tleroy@bird2login.univ-nantes.fr:~/Transferts/TSO
/sandbox/users/tleroy/Transferts/POP_CDS
