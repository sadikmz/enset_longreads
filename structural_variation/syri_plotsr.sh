#!/bin/bash

ref=wildb
qry=wildc
ref_genome=../../${ref}.fna
qry_genome=../../${qry}.fna
EG=../../Ensete_glaucum.fna
out_prefix=${qry}_${ref}
threads=48

minimap2 -ax asm5 -t $threads -k 19 --eqx $ref_genome $qry_genome | samtools view - -Sb -@$threads | samtools sort - -@$threads -o {out_prefix}.bam

ref=wildc
qry=epo
minimap2_out=${qry}_${ref}
minimap2 -ax asm5 -t $threads -k 19 --eqx $ref_genome $qry_genome | samtools view - -Sb -@$threads | samtools sort - -@$threads -o {out_prefix}.bam

ref=epo
qry=mazia
minimap2_out=${qry}_${ref}
minimap2 -ax asm5 -t $threads -k 19 --eqx $ref_genome $qry_genome | samtools view - -Sb -@$threads | samtools sort - -@$threads -o {out_prefix}.bam

ref=mazia
qry=EG
minimap2_out=${qry}_${ref}
minimap2 -ax asm5 -t $threads -k 19 --eqx $ref_genome $EG | samtools view - -Sb -@$threads | samtools sort - -@$threads -o {out_prefix}.bam

# Running syri for finding structural rearrangements between enset assemblies and, EV and EG

syri -c wildc_wildb.bam -r ../../wildb.fna -q ../../wildc.fna -F B --prefix wildc_wildb_
syri -c epo_wildc.bam -q ../../wildc.fna -r ../../epo.fna -F B --prefix epo_wildc_
syri -c mazia_epo.bam -r ../../epo.fna -q ../../mazia.fna -F B --prefix mazia_epo_
syri -c EG_mazia.bam -q ../../mazia.fna -r ../../EG.fna -F B --prefix EG_mazia_

# filter out SNP/HDR/TDM/INS

grep -v -E "CPG|CPL|DEL|HDR|INS|SNP|TDM" wildc_wildb_syri.out > wildc_wildb_syri.filtered.out 
grep -v -E "CPG|CPL|DEL|HDR|INS|SNP|TDM" epo_wildc_syri.out > epo_wildc_syri.filtered.out 
grep -v -E "CPG|CPL|DEL|HDR|INS|SNP|TDM" mazia_epo_syri.out > mazia_epo_syri.filtered.out 
grep -v -E "CPG|CPL|DEL|HDR|INS|SNP|TDM" EG_mazia_syri.out > EG_mazia_syri.filtered.out 

# Run plotsr 
plotsr \
--sr wildc_wildb_syri.filtered.out   \
--sr epo_wildc_syri.filtered.out    \
--sr mazia_epo_syri.filtered.out  \
--sr EG_mazia_syri.filtered.out    \
--genomes genome_file_EV_EG.txt \
--cfg base.cfg \
--itx \
-R \
-H 4 \
-W 7 \
-d 1000 -o EV_EG_syri_minimap2_plotsr.pdf


plotsr \
--sr wildc_wildb_syri.filtered.out   \
--sr epo_wildc_syri.filtered.out    \
--sr mazia_epo_syri.filtered.out  \
--genomes genome_file_EV.txt \
--cfg base.cfg \
--itx \
-R \
-H 4 \
-W 7 \
-d 1000 -o EV_syri_minimap2_plotsr.pdf
