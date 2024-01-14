#!/bin/bash

module purge
#module load intel impi imkl
module load GCC/9.3.0
module load OpenMPI/4.0.3


genotype= #mazia/wildb/wildc/epo
mkdir ${genotype}_output
outdir=${genotype}_output
cd $outdir

ref_genome=${genotype}.fna  
pacbio_hifi=/data/${genotype}.mito_chlr.excluded.adapt.discarded.fasta.gz
PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
SAMTOOOLS=/home/lifesci/lfrwtp/apps/samtools-1.18
out_name=${genome}_ngmlr
threads=48

# run ngmlr 
ngmlr \
-t $threads \
-r $ref_genome \
-q $pacbio_hifi \
-x pacbio \
-o ${out_name}.sam

# convert sam to bam and sort bam 
$SAMTOOOLS/samtools view -Sb -@$threads ${out_name}.sam | $SAMTOOOLS/samtools sort  -@$threads -o ${out_name}.sorted.bam

# index 
$SAMTOOOLS/samtools index ${out_name}.sorted.bam 


# run sniffles

sniffles \
--input ${out_name}.sorted.bam \
--reference ${ref_genome} \
--threads ${threads} \
--minsupport 100 \
--minsvlen 60 \
--vcf ${out_name}.vcf \
--sample-id mz \
--output-rnames  

cat mazia_output/mazia_ngmlr.vcf | sniffles2_vcf_parser.py parsesv > ngmlr_sniffle_stvar_mazia.txt
cat epo_output/epo_ngmlr.picard.vcf | sniffles2_vcf_parser.py parsesv > ngmlr_sniffle_stvar_epo.txt
cat wildb_output/wildb_ngmlr.picard.vcf | sniffles2_vcf_parser.py parsesv > ngmlr_sniffle_stvar_wildb.txt
cat wildc_output/wildc_ngmlr.picard.vcf | sniffles2_vcf_parser.py parsesv > ngmlr_sniffle_stvar_wilc.txt
