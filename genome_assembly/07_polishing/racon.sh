#!/bin/bash
#SBATCH --job-name=polish
#SBATCH --partition=cnode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module purge
module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281

#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREADS

map=map-pb
cpus=28
reads=$(cat input.fofn)
genotype= 
out=${genotype}.1
ref=${genotype}.fna



# Collect repetitive 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

# Map and sort reads
winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t $cpus $ref $reads > $out.sam
samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam

# Merge and index
samtools merge -O bam -@$cpus ${genotype}.bam *.sort.bam
samtools index ${genotype}.bam

# Filter to get primary read alignments
samtools view -F0x104 -@$cpus -hb ${genotype}.bam > ${genotype}.pri.bam
samtools index ${genotype}.pri.bam

## extract coverage
~/apps/T2T-Polish/coverage/sam2paf.sh ${genotype}.pri.bam ${genotype}.pri.paf ${genotype}_coverage1x

racon -t 28 \
$reads  \
${genotype}.pri.paf   \
$ref > $out.polished.fna

rm -rf merylDB repetitive_k15.tx $out.sam 

## round2
ref=${genotype}.1.polished.fna
map=map-pb
cpus=28
reads=$(cat input.fofn)
out=${genotype}.2
prefix=epo_psedochr.2

# Collect repetitive 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

# Map and sort reads
winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t $cpus $ref $reads > $out.sam
samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam

# Merge and index
samtools merge -O bam -@$cpus ${genotype}.bam $out.sort.bam
samtools index ${genotype}.bam

# Filter to get primary read alignments
samtools view -F0x104 -@$cpus -hb ${genotype}.bam > ${genotype}.pri.bam
samtools index ${genotype}.pri.bam

## extract coverage
~/apps/T2T-Polish/coverage/sam2paf.sh ${genotype}.pri.bam ${genotype}.pri.paf ${genotype}_coverage_racon2x

## polish with racon
racon -t 28 \
$reads  \
${genotype}.pri.paf   \
$ref > $out.polished.fna

rm -rf merylDB repetitive_k15.tx $out.sam 

##round3

ref=${genotype}.2.polished.fna
map=map-pb
cpus=28
reads=$(cat input.fofn)
out=${genotype}.3
prefix=epo_psedochr.3

# Collect repetitive 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

# Map and sort reads
winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t $cpus $ref $reads > $out.sam
samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam

# Merge and index
samtools merge -O bam -@$cpus ${genotype}.bam $out.sort.bam
samtools index ${genotype}.bam

# Filter to get primary read alignments
samtools view -F0x104 -@$cpus -hb ${genotype}.bam > ${genotype}.pri.bam
samtools index ${genotype}.pri.bam

## extract coverage
~/apps/T2T-Polish/coverage/sam2paf.sh ${genotype}.pri.bam ${genotype}.pri.paf ${genotype}_coverage_racon2x

## Polish with racon
racon -t 28 \
$reads  \
${genotype}.pri.paf   \
$ref > $out.polished.fna
