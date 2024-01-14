#!/bin/bash

module purge
module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281
MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$MY_NUM_THREADS


# convert bam subreads to fasta

bam2fasta data/subreads.bam -o EV_epo.subreads.pacbio_clr

pacbio_clr=data/EV_epo.subreads.pacbio_clr.fasta
# canu 

canu -assemble -haplotype \
-p EV_wild_epo -d EV_wild_epo \
genomeSize=570m \
UseGrid=true \
gridOptions="--time=48:00:00" \
corMinCoverage=4 \
corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
corMhapFilterThreshold=0.0000000002 \
minOverlapLength=100 \
corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" \
mhapMemory=60g \
mhapBlockSize=500 \
ovlMerDistinct=0.975 \
utgOvlErrorRate=0.065 \
correctedErrorRate=0.045 \
corErrorRate=0.25 \
enableOEA="true" \
-pacbio $pacbio_clr 2>EV_wild_epo.asm.log