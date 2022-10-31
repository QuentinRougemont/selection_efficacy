
#verfiy that all XP start with ATG:
for i in fasta_files_withoutquantiles_CDS/XP_*.fst ; do 
    sed -n 2p $i|\
    sed 's/./& /g' |\
    cut -f 1-3 -d " " >> verify_atg;
done

for i in fasta_files_withoutquantiles_CDS/XP_*.fst ; do 
    sed -n 2p $i|\
    sed 's/./& /g' |\
    awk '{print $(NF-2),"\t",$(NF-1),"\t",$NF}' >> verify_STOP_codon;
done  


mkdir LENGTH
for i in fasta_files_withoutquantiles_CDS/XP*fst ; do ./awk_fasta_length.sh $i ; done
for i in  LENGTH/*txt ; do tail -n 1 $i |awk -v var=$i '{print var"\t"$1 }'>> comptage ; done 



#probleme avec cleanAlignement, nécessaire recompilé!
https://stackoverflow.com/questions/6941332/anticipate-kernel-too-old-errors-between-2-6-16-and-2-6-26-kernel-versions

mkdir 00_data

mkdir fasta_files_withoutquantiles_scaffold

for i in 01_scripts/02_vcf2fasta_NW* ; do msub $i  ; done 
for i in 01_scripts/02_vcf2fasta_NC_0342* ; do msub $i  ; done 
for i in 01_scripts/02_vcf2fasta_NC_03419* ; do msub $i  ; done 


awk '$3=="CDS" {print $9}' gfffile|cut -d ";" -f1 |sort|uniq |wc -l


awk '$3=="CDS" {print $0 }' gfffile |cut -d ";" -f 1 |cut -f 1,4,5,7-9 |sort -u -k6 |sed 's/ID=cds-//g' > sorted_list_cds_uniq.txt



library(dplyr)
#/project/lbernatchez/users/qurou/09.epic4_onV2FULL/01.analysis_refait/GATK/16-wgs_filter/PAL

a=read.table("sorted_successfullXP") #list of XP that were outdputted to CDS folder
b=read.table("sorted_list_cds_uniq.txt")  #all CDS 
colnames(b)[6] <- "XP"
w=anti_join(b,a) 

NC92 <- subset(w, w$V1=="NC_034192.2")     
NC92 <- NC92[order(NC92$V2), ] 

 NC92$length=NC92$V3-NC92$V2  
 
 
## travail dans /project/lbernatchez/users/qurou/09.epic4_onV2FULL/methode_thibaults/BER/NC92:
sed -n 2p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind1.1.txt
sed -n 4p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind1.2.txt
sed -n 6p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind2.1.txt
sed -n 8p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind2.2.txt
sed -n 10p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind3.1.txt
sed -n 12p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind3.2.txt
sed -n 14p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind4.1.txt
sed -n 16p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind4.2.txt
sed -n 18p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind5.1.txt
sed -n 20p fasta_files_withoutquantiles_scaffold/NC_034192.2.fst |grep -o . >> all_position.ind5.2.txt
paste all_position.ind* >> position.tmp
seq 68265699 > snpseq 
paste snpseq position.tmp > position_sequence.txt

awk '$3=="CDS" {print $0 }' gff.92 |cut -d ";" -f 1 |cut -f 1,4,5,7-9 >> CDS.info.txt

