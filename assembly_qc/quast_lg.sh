#!/bin/bash
#SBATCH --job-name=wcquast
#SBATCH --partition=cnode
#SBATCH --nodes=16
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5012
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk
module load intel impi imkl

# MY_NUM_THREADS=$SLURM_CPUS_PER_TASK

# export OMP_NUM_THREADS=$MY_NUM_THREADS

genotype= #mazia/wildb/wildc/epo

pacbio_hifi=absolute_path_to_${genotype}.mito_chlr.excluded.adapt.discarded.fasta.gz
pacbio_clr=absolute_path_to_${genotype}.canu_trimmed.fasta.gz # for Epo

quast-lg.py \
data_${genotype}/${genotype}_primary/${genotype}_s33_k61_adapt_discarded_hic.fasta \
data_${genotype}/${genotype}_pseduochr/EV_${genotype}.final.fna \
data_${genotype}/${genotype}_racon1/${genotype}_out.1.polished.fna \
data_${genotype}/${genotype}_racon2/${genotype}_out.2.polished.fna \
data_${genotype}/${genotype}_racon3/${genotype}_out.3.polished.fna \
data_${genotype}/${genotype}_racon3_pilon1/EV_${genotype}_v3.softmasked.fna \
-o quastouputt \
--split-scaffolds \
-L \
--eukaryote \
--large \
--extensive-mis-size 15000 \
--min-contig 1000 \
--min-alignment 500 \
--k-mer-stats \
--k-mer-size 17 \
--circos \
--use-all-alignments \
--plots-format png \
--threads 28 \
--pacbio $pacbio_hifi \
--report-all-metrics







