#!/bin/bash

# unknown\unspecified repeatmasker output 
unknown_unspecified.RM.out=${genotype}.unknown_unspecified.RM.out.fna
#repeatmasker final output
rmout=${genotype}.RM.out.fasta


## extract unknown or unspecified repeats
rmsk2bed < ${genome}.fna.out | cut -f1-3,11 | sort -V | uniq > ${genome}.fna.out 
bedtools getfasta -fi ${genome}.fna -bed ${genome}.fna.out.bed > ${genome}.fna.out.fasta 

# annotate repetitve sequences identified by RepeatMasker using a custome repeat library ${genotype}.repeat_lib.fa generated in repeat_analysis.sh

BASENAME=$(basename $rmout | sed 's/.RM.final.out.fasta//g')

dante \
-q ${BASENAME}.RM.final.out.fasta \
-D Viridiplantae_v3.0 \
-o dante_${BASENAME}_rmout.gff \
-M BL62 \
-c 64

dante_gff_output_filtering.py -dg dante_${BASENAME}_rmout.gff -dir dante_${BASENAME}_out


