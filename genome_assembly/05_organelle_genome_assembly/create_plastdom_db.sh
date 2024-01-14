#!/bin/bash

# To create the custom organellar database, get_organelle assembled enset mitochondrial sequences were concatenated with reference mitochondrial sequences from NCBI (https://ftp.ncbi.nlm.nih.gov/refseq/release/ mitochondrion/), mitochondrial sequences of Phoenix dactylifera (NC_016740.1), and Musa acuminata. 

path_data=/home/primary_assembly/
genotype= #

# makeblastdb -dbtype nucl -in mito.db/mitoRefSeq_MA_PD.mito.fasta 

makeblastdb -dbtype nucl -in ${path_data}/${genotype}.fna 

blastn \
-query ${path_data}/${genotype}.fna  \
-db mito.db/mitoRefSeq_MA_PD.mito.fasta  \
-outfmt 6 'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out ${genotype}.refmito.e99.blastn.out 

Rscript parse_blastn.R


blastn \
-db ${path_data}/${genotype}.fna  \
-query mito.db/mitoRefSeq_MA_PD.mito.fasta  \
-outfmt 6 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out refmito.${genotype}.e99.blastn.out 


Rscript parse_blastn.R

# Downloaded mitochorial sequences maching >80 % coverage and identify with sequences in EV genomes  

python download_mitoseq_genbank.py -i rifmito_ID_EV.txt -o rifmito_ID_EV.80_80.gb

# convert genbank into fasta and concatenate all mito-sequences into one fasta file
cat rifmito_ID_EV.80_80.fasta NC_016740.fast MA.mito.fasta EV_get_organelle.fasta > mito.db/mitoRefSeq_MA_PD_EV_chrpt.mitofasta