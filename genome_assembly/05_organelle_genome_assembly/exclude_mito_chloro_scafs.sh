#!/bin/bash

genotype= #mazia/wildb/wildc/epo
asm_data=/home/data/primary_assembly/${genotype}.fna
out_dir=blast_out

makeblastdb -dbtype nucl -in mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa
makeblastdb -dbtype nucl -in ${asm_data} 

blastn \
-query ${asm_data} \
-db mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa \
-outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out blast_out/${genotype}.mito_chl.e99.blastn.out 

blastn \
-db ${asm_data} \
-query mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa \
-outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out blast_out/refmito.${genotype}.e99.blastn.out

# swap bed coordinates that start position > end
awk '{if ($9 > $10) print $1"\t"$10"\t"$9; else print $1"\t"$9"\t"$10}' blast_out/${genotype}.mito_chl.e99.blastn.out > blast_out/${genotype}.q.blastn.bed
awk '{if ($11 > $12) print $2"\t"$12"\t"$11; else print $2"\t"$11"\t"$12}' blast_out/refmito.${genotype}.e99.blastn.out > blast_out/${genotype}.r.blastn.bed

#awk '{if ($9 > $10) print $1"\t"$10"\t"$9"\t"$3; else print $1"\t"$9"\t"$10"\t"$3}' blast_out/${i}.e99.blastn.out > blast_out/${i}.q.blastn.bed

# combine 
cat blast_out/${genotype}.q.blastn.bed blast_out/${genotype}.r.blastn.bed > blast_out/${genotype}.blastn.bed

# sort bed files
bedtools sort -i blast_out/${genotype}.blastn.bed > blast_out/${genotype}.blastn.sorted.bed

# Merge coordinated
bedtools merge -i blast_out/${genotype}.blastn.sorted.bed > blast_out/${genotype}.blastn.sorted.merged.bed

# #
#samtools faidx ${asm_data}
cat ${asm_data}.fai | awk '{print $1, $2}' OFS='\t' > blast_out/${genotype}.contig.sizes

Rscript mito_chlo_contig_cov_80.R

grep -v -f blast_out/${genotype}_mito_chrpt_contigs.cov_80.txt blast_out/${genotype}.contig.sizes | awk '{print $1}' > blast_out/${genotype}_non_plastome_contigs.txt 
seqtk subseq ${asm_data} blast_out/${genotype}_non_plastome_contigs.txt > blast_out/${genotype}.fna 