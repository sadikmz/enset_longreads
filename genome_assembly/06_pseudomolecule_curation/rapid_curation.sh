#!/bin/bash

# Manual curation of Omni-C scaffolds was performed using rapid curation pipeline (https://gitlab.com/wtsi-grit/rapid-curation)

genotype=wildb
primary_assembly=/home/data/assembly/assembled_fasta/${genotype}_s33_k61.fasta
ln -s $primary_assembly .

# source https://www.youtube.com/watch?v=LWy6pwCQNDU

bash ~/apps/rapid-curation/sif/rc_suite_orac.sh -p hic,coverage,telomere,gap,repeat -s TTTAGGG -t mazia

## Find telomeric sequences

#option #1

genotype=ref

tidk search ${genotype}.fa -s TTTAGGG -o ${genotype}.T3AG3.tidk --dir tidk_out --extension tsv 
tidk search ${genotype}.fa -s CCCTAAA -o ${genotype}.C3TA3.tidk --dir tidk_out --extension tsv 


# Option #2
samtools faidx ${genotype}.fa
cat ${genotype}.fa.fai | awk '{print $1, $2}' OFS='\t' > ${genotype}.genome_size.txt
bedtools makewindows -g ${genotype}.genome_size.txt -w 10000 > ${genotype}_windows.bed 

bedtools nuc -fi ${genotype}.fa -bed ${genotype}_windows.bed -pattern TTTAGGG -seq | grep -v '#' | cut -f1,2,3,13 > ${genotype}.10k.T3AG3.bed
bedtools nuc -fi ${genotype}.fa -bed ${genotype}_windows.bed -pattern CCCTAAA -seq | grep -v '#' | cut -f1,2,3,13 > ${genotype}.10k.C3TA3.bed

## embed coverage, repeat ... into scaffold-hic interaction map 

zcat bedgraph.file.gz | PretextGraph -i input.pretext -n "graph name"
bigWigToBedGraph bigwig.file /dev/stdout | PretextGraph -i input.pretext -n "graph name"

~/apps/rapid-curation/rapid_split.pl -fa ${genotype}_s33_k61.fasta

# Vizualize scaffold-hic interaction contact map and anchor, order, orient and correct misassemblies of scaffolds into pseudochromosomes using the genomic proximity signal of the Omni-C reads

~/apps/rapid-curation/rapid_pretext2tp_XL.py scaffold_divided.tpf original.pretext.agp_1 
~/apps/rapid-curation/rapid_join.pl -fa ref.fa -tpf rapid_prtxt_XL.tpf -csv chrs.csv -out ${genotype}_rc -hap haps_rapid_prtxt_XL.tpf 







