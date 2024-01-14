## Pilon 

genotype= 
Fread=/home/data/EV_NGS/DP8400009105TR_L01_520_1_val_1.fq.gz
Rread=/home/data/EV_NGS/DP8400009105TR_L01_520_2_val_2.fq.gz
ref=${genotype}_out.3.polished.fna # racon 3x polished 
cpus=28
out=${genotype}_p1
bam=${out}.sorted.bam 


# round1
bwa-mem2 index $ref
bwa-mem2 mem -t 28 $ref $Fread $Rread | samtools view - -Sb -@28 | samtools sort - -@28 -o $out.sorted.bam
pilon --genome $ref --frags $bam --fix all --output pilon1 --outdir round1 --changes --diploid --tracks --vcfqe --verbose

# Pilon round
out=round1/${genotype}_p2
ref=pilon1.fasta 
bwa-mem2 mem -t 28 $ref $Fread_CL $Rread_CL | samtools view - -Sb -@28 | samtools sort - -@28 -o $out.sorted.bam
samtools merge -O bam -@28 --write-index -o ${genotype}_merged.bam ${genotype}_p1.sorted.bam ${genotype}_p2.sorted.bam 

#run pilon 
pilon --genome $ref --frags $bam --fix all --output pilon2 --outdir round2 --changes --diploid --tracks --vcfqe --verbose


# Pilon round3 
out=round1/${genotype}_p3
ref=round2/pilon2.fasta 
bwa-mem2 mem -t 28 $ref $Fread_CL $Rread_CL | samtools view - -Sb -@28 | samtools sort - -@28 -o $out.sorted.bam
samtools merge -O bam -@28 --write-index -o ${genotype}_merged.bam ${genotype}_p1.sorted.bam ${genotype}_p2.sorted.bam 

#run pilon 
pilon --genome $ref --frags $bam --fix all --output pilon3 --outdir round3 --changes --diploid --tracks --vcfqe --verbose
