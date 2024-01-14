### Epo

genotype=epo
epo_pacbio_clr=canu_corrected_${genotype}_correctedReads.fasta.gz
threads=48
genome_size=550m
output=epo_clr


# Flye 
flye \
--pacbio-corr $epo_pacbio_clr \
--out-dir $output \
--iterations 3 \
--keep-haplotypes \
--genome-size $genome_size \
--min-overlap 8000 \
--scaffold \
--threads $threads




