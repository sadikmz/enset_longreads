#!/bin/bash
#SBATCH --job-name=mzmito
#SBATCH --partition=cnode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

module purge
module load intel/2019.5.281-GCC-8.3.0-2.32 impi/2018.5.288 imkl/2019.5.281

MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$MY_NUM_THREADS


#asm_data=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/mazia_hic/assembled_fasta/mazia_s33_k61_adapt_discarded_hic.fasta
asm_data=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/mazia_hic/assembled_fasta/mazia_s33_k61_adapt_discarded_hic.fasta

genome=mazia
out_dir=blast_out

#makeblastdb -dbtype nucl -in mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa
makeblastdb -dbtype nucl -in ${asm_data} 

blastn \
-query ${asm_data} \
-db mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa \
-outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out blast_out/${genome}.mito_chl.e99.blastn.out 

blastn \
-db ${asm_data} \
-query mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa \
-outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
-evalue 1e-99 \
-max_hsps 1 \
-num_threads 28 \
-out blast_out/refmito.${genome}.e99.blastn.out

# swap bed coordinates that start position > end
awk '{if ($9 > $10) print $1"\t"$10"\t"$9; else print $1"\t"$9"\t"$10}' blast_out/${genome}.mito_chl.e99.blastn.out > blast_out/${genome}.q.blastn.bed
awk '{if ($11 > $12) print $2"\t"$12"\t"$11; else print $2"\t"$11"\t"$12}' blast_out/refmito.${genome}.e99.blastn.out > blast_out/${genome}.r.blastn.bed

awk '{if ($9 > $10) print $1"\t"$10"\t"$9"\t"$3; else print $1"\t"$9"\t"$10"\t"$3}' blast_out/${i}.e99.blastn.out > blast_out/${i}.q.blastn.bed

# combine 
cat blast_out/${genome}.q.blastn.bed blast_out/${genome}.r.blastn.bed > blast_out/${genome}.blastn.bed

# sort bed files
bedtools sort -i blast_out/${genome}.blastn.bed > blast_out/${genome}.blastn.sorted.bed

# Merge coordinated
bedtools merge -i blast_out/${genome}.blastn.sorted.bed > blast_out/${genome}.blastn.sorted.merged.bed

# #
#samtools faidx ${asm_data}

cat ${asm_data}.fai | awk '{print $1, $2}' OFS='\t' > blast_out/${genome}.contig.sizes

Rscript mito_chlo_contig_cov_80.R

grep -v -f blast_out/${genome}_mito_chrpt_contigs.cov_80.txt blast_out/${genome}.contig.sizes > blast_out/${genome}_non_plastome_contigs.txt 
seqtk ${asm_data} blast_out/${genome}_non_plastome_contigs.txt > blast_out/${genome}.adapt_plastome_discarded.fasta 


asm_data=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/mazia_hic/assembled_fasta/mazia_s33_k61_adapt_discarded_hic.fasta
genome=mazia

grep -v -f blast_out/${genome}_mito_chrpt_contigs.cov_80.txt blast_out/${genome}.contig.sizes > blast_out/${genome}_non_plastome_contigs.txt 
seqtk ${asm_data} blast_out/${genome}_non_plastome_contigs.txt > blast_out/${genome}.adapt_plastome_discarded.fasta 







