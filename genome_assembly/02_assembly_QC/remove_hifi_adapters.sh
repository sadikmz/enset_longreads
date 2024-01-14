#!/bin/bash

# extract hifi bam subreads 
extracthifi raw.bam hifi.bam 

# convert hifi_bam subreads to fasta 

hifi_bam=absolute_path_to_bam.hifi.subreads.bam
hifi_fasta=absolute_path_to_bam_hifi_raw 
genotype= #mazia/wildb/wildc/ 

bam2fasta $hifi_bam -o $hifi_fasta -j 48 

# adapter removal

cutadapt \
-b "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT;min_overlap=35" \
-b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=35" \
--discard-trimmed \
--revcomp \
-e 0.1 \
--report=minimal \
--discard  \
--output ${genotype}.adapt.discarded.fasta.gz  \
--cores=32 \
${hifi_fasta}.fasta.gz



# mkdir jellyfish_histo

for k in 17 21 27 31 37 41; 
do 
  # This will generate kmer count and kmer histogram for esimated genome size of 580M and with estimated 43X coverage.
  # The  estimated to be gerated kmer size is to accomodate Bedadeti 42?x coverage reads as all other are < 30X 
  do
  ${prog}/jellyfish-linux count -C -m $k -s 570M --bf-size 50G -t 24 <(zcat ${genotype}_hifi.adapt.discarded.fasta.gz) -o ${k}_mer.reads.jf
  ${prog}/jellyfish-linux histo -t 24 ${k}_mer.reads.jf > jellyfish_histo/${k}_mer.histo
  rm -rf ${k}_mer.reads.jf

done




