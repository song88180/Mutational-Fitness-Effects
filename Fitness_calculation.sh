#!/bin/bash

#Mapping reads to the WT reference sequence(150bp target sequence with 6-bp wildtype sequence at each of two ends). I took ADA2 analysis as an example and the rest 20 genes were analyzed in the same way.
bwa mem -t 4 ./ADA2.fasta ./ADA2-t0_R1_001.fastq.gz ./ADA2-t0_R2_001.fastq.gz > ADA2_t0.sam

#Read1 and read2 encode the same mutations and they are fully mapped to the 162bp reference sequence.
awk -F ' ' '{print $1, $9, $13}' *.sam | sed 's/-//g' | sort | uniq -c | awk '$1 > 1 { print $2, $3, $4 }' > 01.paired.same.txt
awk '$2 == "162" {print $0}' 01.paired.same.txt > 02.paired.same.txt

#Excluding indel mutations.
grep '\^' 02.paired.same.txt > 03.indel.txt
comm -1 -3 03.indel.txt 02.paired.same.txt > 04.mutations.txt

#Extracting the reads which are WT and single mutants
awk '{print $9, $10}' FPAT='[0-9]+' 04.mutations.txt | awk  '{ print $1 + $2}' > 05.sumofmatchedbases.txt
paste 04.mutations.txt 05.sumofmatchedbases.txt | awk '$4 >= 161' > 06.onlytheWTandsinglemutations.txt
awk '{print $1}' 06.onlytheWTandsinglemutations.txt | sed 's/NGSNJ086/NGSNJ-086/g' >  07.readswith0or1mutation.txt
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' 07.readswith0or1mutation.txt *.sam > 08.theselectedsamreads.sam

#IIdentifying mutation position, WT base and mutation base
awk '{print $1, $10, $13}' 08.theselectedsamreads.sam | awk 'NR % 2 == 1' | grep 'TGCTGACCATATAGGCAGCAGAGGCAAAGAAGAAGTTAAG' | grep 'CCCATGGCATC' > 09.thesequenceswithintacsurrongdingsequences.txt
awk '{print $2}' 09.thesequenceswithintacsurrongdingsequences.txt | awk -F 'GGCAAAGAAGAA' '{print $2}' | awk '{print substr($1, 1, 162)}' > 091.the_sequencesof162segments.txt
paste -d' ' 09.thesequenceswithintacsurrongdingsequences.txtt 091.the_sequencesof162segments.txt | awk '{print $1, $4, $3}' | awk '{print $8}' FPAT='[0-9]+' | awk '{print $1 = $1 +1}' > 092.thepositionofmutation.txt
paste -d' ' 09.thesequenceswithintacsurrongdingsequences.txtt 091.the_sequencesof162segments.txt | awk '{print $1, $4, $3}' | awk '{print $6}' FPAT='[A-Z]+' > 093.printtheWTbase.txt
paste -d' ' 09.thesequenceswithintacsurrongdingsequences.txtt 091.the_sequencesof162segments.txt | awk '{print $1, $4, $3}' > 094.printthenamesequenceandMD.txt
paste -d' ' 094.printthenamesequenceandMD.txt 092.thepositionofmutation.txt 093.printtheWTbase.txt  > 095.printthemutationpositionandWTbaseinafile.txt
awk '{print substr($2, $4, 1)}' 095.printthemutationpositionandWTbaseinafile.txt  > 096.themutatedbase.txt

#Counting the WT genotype and single-mutation genotypes which are covered by more than 50 reads
paste -d' ' 095.printthemutationpositionandWTbaseinafile.txt 096.themutatedbase.txt | awk '{print $0, $4"-"$5"-"$6}' | awk '{print $2, $7}'| sort -k2,2 | uniq -c | awk '{if ( $1 >= 50 ) print }' | sed -e '1s/$/WT/' | grep 'CCCATG' > 097.allthegenotypesandthenumber.txt

#Calculating the ratios of mutant genotype divided by WT genotype
awk 'BEGIN{first_line=0;divide_by=1;}{if(first_line==0){first_line=1; divide_by=$1;}print $1/divide_by;}END{}' 097.allthegenotypesandthenumber.txt | paste - 097.allthegenotypesandthenumber.txt > 098.allthegenotypesandtheratio.txt


#Calculating genotype fitness based on the genotype frequency change between t0 and a later time point t. N is the WT generation number between two time points. 
join -1 1 -2 1 <(sort -k 1 098.allthegenotypesandtheratio.txt) <(sort -k 1 Path_to_another_time_point/098.allthegenotypesandtheratio.txt) | awk '{print $1, $2, $5}' | awk '{print $1, ($3/$2)^(1/N)}’



