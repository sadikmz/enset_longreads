#!/bin/bash


# Installation 

conda install -c bioconda extracthifi pbindex 

conda update extracthifi

# a2_EV_Mazia
# a3_EV_Wild_B
# a4_EV_Wild_C

# extract hifi bam subreads : extracthifi is used to extract PacBio HiFi reads (>= Q20) from full CCS output (reads.bam).
# extracthifi raw.bam hifi.bam 

extracthifi X204SC21052271-Z01-F001_03/raw_data/a4/m64164_210626_142631.bc1068.subreads.bam  X204SC21052271-Z01-F001_03/raw_data/a4/m64164_210626_142631.bc1068.hifi.subreads.bam

# extract reads that did not make make it Q>=20
bamtools filter \
-in X204SC21052271-Z01-F001_03/raw_data/a4/m64164_210626_142631.bc1068.subreads.bam \
-out X204SC21052271-Z01-F001_03/raw_data/a4/a4_hifi_reads_less_than_Q20.bam -tag "rq":"<0.99"

for i in a2 a3
bamtools filter \
-in X204SC21051840-Z01-F001/raw_data/${i}/m64164_210627_204210.subreads.bc1015--bc1015.bam \
-out X204SC21051840-Z01-F001/raw_data/${i}/${i}_hifi_reads_less_than_Q20.bam -tag "rq":"<0.99"

# https://github.com/marbl/merqury/issues/76
# samtools view reads.ccs.bam |awk '{RQ=-1; for (i = NF; i > 11; i--) { if (match($i, "rq:")) RQ=substr($i, 6, length($i)); } if (RQ >= 0.99) { print "@"$1; print $10; print "+"; print $11; }}'  > reads.ccs.fastq
# # or 
# samtools view reads.ccs.bam  |awk '{RQ=-1; for (i = NF; i > 11; i--) { if (match($i, "rq:")) RQ=substr($i, 6, length($i)); } if (RQ >= 0.99) { print ">"$1; print $10; }}' |fold -c > reads.ccs.fasta


bam2fasta X204SC21051840-Z01-F001/raw_data/${i}/${i}_hifi_reads_less_than_Q20.bam -o X204SC21051840-Z01-F001/raw_data/${i}/${i}_hifi_reads_less_than_Q20


# on avon 
bamtools filter \
-in /home/lifesci/lfrwtp/data/PacBio_Hifi/X204SC21051840-Z01-F001/raw_data/a2/m64164_210627_204210.subreads.bc1015--bc1015.bam \
-out /home/lifesci/lfrwtp/data/PacBio_Hifi/X204SC21051840-Z01-F001/raw_data/a2/m64164_210627_204210.subreads.bc1015--bc1015_less_Q20.bam -tag "rq":"<0.99"

bam2fasta /home/lifesci/lfrwtp/data/PacBio_Hifi/X204SC21051840-Z01-F001/raw_data/a2/m64164_210627_204210.subreads.bc1015--bc1015_less_Q20.bam \
 -o /home/lifesci/lfrwtp/data/PacBio_Hifi/X204SC21051840-Z01-F001/raw_data/a2/mazia_less_Q20


done 

cutadapt \
-b "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT;min_overlap=35" \
-b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=35" \
--discard-trimmed \
--revcomp \
-e 0.1 \
--report=minimal \
--discard  \
--output wild_c_hifi.adapt.discarded.fasta.gz  \
--cores=32 \
wild_c_hifi.fasta.gz



# convert hifi_bam subreads to fasta 

bam2fasta /home/u1866313/mnt/Shared269/MGrant/PacBio_HiFi/X204SC21052271-Z01-F001_multipath/X204SC21052271-Z01-F001_01/HiFi_data/a4/m64164_210626_142631.bc1068.hifi.subreads.bam -o /home/u1866313/mnt/Shared269/MGrant/PacBio_HiFi/X204SC21052271-Z01-F001_multipath/X204SC21052271-Z01-F001_01/HiFi_data/a4/wild_c_hifi

# Generate Kmer histogram
prog=/home/u1866313/apps/jellyfish/

ln -s /home/u1866313/mnt/Shared269/MGrant/PacBio_HiFi/X204SC21052271-Z01-F001_multipath/X204SC21052271-Z01-F001_01/HiFi_data/a4/wild_c_hifi.fasta.gz .


# adapter removal

cutadapt \
-b "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT;min_overlap=35" \
-b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=35" \
--discard-trimmed \
--revcomp \
-e 0.1 \
--report=minimal \
--discard  \
--output wild_c_hifi.adapt.discarded.fasta.gz  \
--cores=32 \
wild_c_hifi.fasta.gz


# QV value 
# build k-mer hash tables for high-coverage reads; discard singletons
./yak count -b37 -t32 -o ccs.yak ccs-reads.fq.gz

# compute k-mer QV for reads
./yak inspect ccs.yak sr.yak > ccs-sr.kqv.txt

# print k-mer histogram
./yak inspect sr.yak > sr.hist


# mkdir jellyfish_histo

for k in 17 21 27 31 37 41; 
do 
  # This will generate kmer count and kmer histogram for esimated genome size of 580M and with estimated 43X coverage.
  # The  estimated to be gerated kmer size is to accomodate Bedadeti 42?x coverage reads as all other are < 30X 
  do
  ${prog}/jellyfish-linux count -C -m $k -s 570M --bf-size 50G -t 24 <(zcat wild_c_hifi.adapt.discarded.fasta.gz) -o ${k}_mer.reads.jf
  ${prog}/jellyfish-linux histo -t 24 ${k}_mer.reads.jf > jellyfish_histo/${k}_mer.histo
  rm -rf ${k}_mer.reads.jf

done




