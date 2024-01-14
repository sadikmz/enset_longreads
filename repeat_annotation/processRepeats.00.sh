#!/bin/bash


ProcessRepeats -species viridiplantae -nolow -noint -a EV_wildb.final.v1.fna.cat.gz -excln -poly 

rmsk2bed --input EV_wildb.final.v1.fna.out --output EV_wildb.final.v1.fna.out.bed --keep-header

rmsk2bed < EV_wildb.final.v1.fna.out --keep-header > EV_wildb.final.v1.fna.out.bed --max-mem 4G 

# TEsorter 

for i in EV_mazia EV_wildb EV_wildc EV_epo 
do
TEsorter \
${i}.unmasked.fna \
--seq-type nucl \
--processors 36 \
--pass2-rule 80-80-80 \
--prefix ${i}_tesorter.out \
-db rexdb-plant \
-genome
done 

path_RMbed=/home/u1866313/hifi_assembly/repeat_annotation/all_RM_out_final

for i in EV_mazia EV_wildb EV_wildc EV_epo 
do

 gff2bed < ${i}_tesorter.out.dom.gff3 > ${i}_tesorter.out.dom.bed
 cat ${i}_tesorter.out.dom.bed | grep -E 'ypsy|opia' > ${i}_tesorter.gypsy_copy.bed
 cat ../all_RM_out_final/${i}.fna.out.bed | cut -f1-3 > ${i}.fna.out.bed
 bedtools intersect -a ${i}.fna.out.bed -b ${i}_tesorter.gypsy_copy.bed -wao |  grep 'TEsorter' > ${i}.rmsk2bed.tesort.intersect.txt; 

done 
