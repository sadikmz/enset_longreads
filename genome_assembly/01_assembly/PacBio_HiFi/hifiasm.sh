#!/bin/bash

# adapter and organelar sequences undiscaded hifi reads 
genotype= # wildb/mazia/wildc
cleaned_hifi_reads=/data/${genotype}_hifi.mito_chlr.excluded.adapt.discarded.fasta.gz
omnic_Reads=/home/data/omnic_Reads

# Assembly PacBio-Hifi reads using Hifiasm 
~/apps/hifiasm-0.17.7/hifiasm  \
-o mazia_s33_k61 \
-t 48 \
--hg-size 570m \
--k 61 \
-s 0.33 \
--h1 ${omnic_Reads}/DTG_OmniC_1.fq.gz \
--h2 ${omnic_Reads}/DTG_OmniC_2.fq.gz  \
$cleaned_hifi_reads


# assembly statistics 

./assembly_stat.sh