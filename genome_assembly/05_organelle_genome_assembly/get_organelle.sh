#!/bin/bash

PATH_WGS=/home/data/EV_NGS_accessions/trim_out


for i in $(cat EV_NGS_accessions.txt)
do
get_organelle_from_reads.py \
-1 PATH_WGS/${i}.1_val_1.fq.gz \
-2 S${i}.2_val_2.fq.gz \
-o ${i}_plastome \
-R 30 \
-F embplant_pt \
--max-n-words 2E9 \
-t 48 \
-k 17,21,55,85,115,127 \
--remove-duplicates 2E7 \
--reverse-lsc \
--reduce-reads-for-coverage inf \
--max-reads inf
done 


for i in $(cat EV_NGS_accessions.txt)
do 

	cat "$i"_plastome/extended_K*.assembly_graph.fastg.extend-embplant_pt-embplant_mt.fastg >> EV_WGS_chloroplast.fasta
	cat "$i"_plastome/extended_K*.assembly_graph.fastg.extend-embplant_pt-embplant_mt.fastg | grep '^>' >> EV_WGS_chloroplast.count.txt

done 