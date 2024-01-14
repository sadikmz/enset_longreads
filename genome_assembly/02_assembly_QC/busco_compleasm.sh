#!/bin/bash

# Primary assembly 
genotype= #mazia/wildb/wildc/epo
prim_asm=/data/prim_asm/${genotype}_s33_k61_adapt_discarded_hic.fasta

# Pseudo_chr_assembly 
pseudo_chr_asm=/data/pseudo_chr_asm/${genotype}.fna 

module purge
module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281

genome=mazia
output_dir=${genotype}_busco_compleasm_out
lineage=embryophyta_odb10
mode=genome
threads=28
dowloaded_lineage_path=compleasm_lineage_dir
mode=

# For primary assemblies 
busco \
-i $pri_asm  \
-m genome  \
-l $lineage \
-c $threads \
-o busco_${genotype}_primasm.out \
--miniprot \
--long \
--tar \
--offline 

compleasm download embryophyta -L compleasm_lineage_dir

compleasm \
run \
-a  $pri_asm \
-l  $lineage \
-t $threads \
-o compleasm_${genotype}_primasm_out \
-L $dowloaded_lineage_path

# For pseudochromosome assemblies 

busco \
-i $pseudo_chr_asm  \
-m genome  \
-l $lineage \
-c $threads \
-o busco_${genotype}_psedochr_out \
--miniprot \
--long \
--tar \
--offline 


compleasm \
run \
-a  $pri_asm \
-l  $lineage \
-t $threads \
-o compleasm_${genotype}_psedochr_out \
-L $dowloaded_lineage_path