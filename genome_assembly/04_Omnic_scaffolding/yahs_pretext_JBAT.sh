#!/bin/bash

## tools setup
## need juicer_tools/pretextmap and samtools if want to do hic plot
## juicer_tools: https://github.com/aidenlab/juicer/wiki/Download
## PretextMap: https://github.com/wtsi-hpag/PretextMap
## PretexSnapshot: https://github.com/wtsi-hpag/PretextSnapshot
## samtools: https://github.com/samtools/samtools
## please adjust the path to juicer_tools and samtools
## here we use 12 CPUs and 32Gb memory for juicer_tools pre - adjust it according to your device
## see more information for juicer tools https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start
## output file will be ${outdir}/${genotype}.hic
## the output hic file could be viewed with JuiceBox https://github.com/aidenlab/Juicebox
#juicer_tools="java -Xmx32G -jar /bin/juicer_tools_1.22.01.jar pre --threads 12"
## v1.9.9 seems much faster than v1.22.01

juicer_tools="java -Xmx32G -jar /home/apps/juicer_out/juicer_tools_1.19.02.jar pre"
pretext_map="/home/miniconda3/envs/hic_pro/bin/PretextMap"
pretext_snapshot="/home/miniconda3/envs/hic_pro/bin/PretextMap/PretextSnapshot"
samtools="/home/miniconda3/envs/hic_pro/bin/samtools"


## input files path and output dir
outdir="out_pretext_JBAT"
genotype="wildb"
scaffold_dir=/home/data/scaffold_dir
contig_dir=/home/data/contigs_dir
omnic_dir=/home/data/omnic_dir
contigs="${contig_dir}/${genotype}.assembly.fna" # need to be indexed, i.e., ${test}.contigs.fasta.gz.fai is presented
hicaln_bycoord="${genotype}.unique_mapped.sorted.bam" # could be .bed, .bam or .bin file
hicaln_bynames="${genotype}.unique_mapped.sorted.bynames.bam" # could be .bed, .bam or .bin file

### generate an alignment 
# source https://omni-c.readthedocs.io/en/latest/fastq_to_bam.html#dups

samtools faidx ${genotype}.fna

bwa-mem2 index -p ${genotype}.fna.index ${genotype}.fna

# align OmniC reads
bwa-mem2 mem -5SP -T0 -t 48 ${genotype}.fna.index ${omnic_dir}/DTG_OmniC_R1.fq.gz  ${omnic_dir}/DTG_OmniC_R2.fq.gz  -o aligned.${genotype}.sam

# Recording valid ligation events
pairtools parse  --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 48 --nproc-out 48 --chroms-path ${genotype}.fna aligned.${genotype}.sam >  parsed.${genotype}.pairsam

# generate alignment stats
samtools flagstat aligned.${genotype}.sam > aligned.${genotype}.sam.all.stats

# remove sam
rm aligned.${genotype}.sam

# Sorting the pairsam file
pairtools sort --nproc 48 --tmpdir tmp/  parsed.${genotype}.pairsam > sorted.${genotype}.pairsam

# remove unsorted pairsam
rm parsed.${genotype}.pairsam 

# Remove PCR duplicates 
pairtools dedup --nproc-in 48 --nproc-out 48 --mark-dups --output-stats stats.txt --output dedup.${genotype}.pairsam sorted.${genotype}.pairsam

# Generating .pairs and bam files
pairtools split --nproc-in 48 --nproc-out 48 --output-pairs mapped.${genotype}.pairs --output-sam unsorted.${genotype}.bam dedup.${genotype}.pairsam

# Generating the final bam file

samtools sort -n -@48 -o mapped.${genotype}.sorted.bynames.bam unsorted.${genotype}.bam

samtools index mapped.${genotype}.sorted.bynames.bam


## run yahs 

# bam sorted by names
yahs \
${contig_dir}/${genotype}.assembly.fna  \
${hicaln_bynames}  \
-l 12000 \
-r 10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000 \
-o ${genotype}_yahs_bynames


## run Juicer 

#### this is to generate input file for juicer_tools - non-assembly mode or for PretextMap
## here we use 8 CPUs and 32Gb memory for sorting - adjust it according to your device
(juicer pre ${scaffold_dir}/${genotype}_yahs_bynames.bin ${scaffold_dir}/${genotype}_yahs_bynames_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T ${outdir} --parallel=28 -S32G | awk 'NF' > ${outdir}/alignments_sorted.txt.part) && (mv ${outdir}/alignments_sorted.txt.part ${outdir}/alignments_sorted.txt)
## prepare chromosome size file from samtools index file
# ${samtools} faidx ${outdir}/${genotype}_scaffolds_final.fa
# cut -f1-2 ${outdir}/${genotype}_scaffolds_final.fa.fai >${outdir}/${genotype}_scaffolds_final.chrom.sizes

## another way to prepare chromosome size file
## this is an easier way especially when we have >2G scaffolds which need scaling 
cat ${outdir}/tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${genotype}_scaffolds_final.chrom.sizes

## do juicer hic map
(${juicer_tools} ${outdir}/alignments_sorted.txt ${outdir}/${genotype}.hic.part ${outdir}/${genotype}_scaffolds_final.chrom.sizes) && (mv ${outdir}/${genotype}.hic.part ${outdir}/${genotype}.hic)

## do Pretext hic map
(awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' ${outdir}/${genotype}_scaffolds_final.chrom.sizes; awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' ${outdir}/alignments_sorted.txt) | ${pretext_map} -o ${outdir}/${genotype}.pretext

# and a pretext snapshot
${pretext_snapshot} -m ${outdir}/${genotype}.pretext --sequences "=full" -o ${outdir}

#### this is to generate input file for juicer_tools - assembly (JBAT) mode (-a)
juicer pre -a -o ${outdir}/${genotype}_JBAT ${scaffold_dir}/${genotype}_yahs_bynames.bin ${scaffold_dir}/${genotype}_yahs_bynames_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre_JBAT.log
cat ${outdir}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${genotype}_JBAT.chrom.sizes
(${juicer_tools} ${outdir}/${genotype}_JBAT.txt ${outdir}/${genotype}_JBAT.hic.part ${outdir}/${genotype}_JBAT.chrom.sizes) && (mv ${outdir}/${genotype}_JBAT.hic.part ${outdir}/${genotype}_JBAT.hic)


#### this is to generate final genome assembly file after manual curation with JuiceBox (JBAT)
## the output assembly file after curation is ${outdir}/${genotype}_JBAT.review.assembly
## the final output is ${outdir}/${genotype}_JBAT.FINAL.agp and ${outdir}/${genotype}_JBAT.FINAL.fa
juicer post -o ${outdir}/${genotype}_JBAT ${outdir}/${genotype}_JBAT.assembly ${outdir}/${genotype}_JBAT.liftover.agp ${contigs}
