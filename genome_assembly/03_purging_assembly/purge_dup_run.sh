#!/bin/bash
#SBATCH --job-name=mz
#SBATCH --partition=devel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=01:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

pri_asm=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/decontamination/mazia.adapt_plastome_discarded_decontam.fasta
pb_list=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/organelles_genome/mazia/minimap2_out/mazia_hifi.mito_chlrpt_adapt_discarded.fasta


genotype= #mazia/wildb/wildc # epo with PacBio-CLR reads
##run busco on assembly 
busco  -i $pri_asm  -m genome  -l embryophyta_odb10 -c 28 -o busco_decontam.out --long -f --offline

# ## extract contigs with duplcated busco genes 

cat busco_decontam.out/run_embryophyta_odb10/full_table.tsv | grep -v "#" | grep 'pt' | grep 'Duplicated' | awk '{print $3}' | sort -V | uniq > ${genotype}_duplicated_busco_contigs.txt

# ## extrac sequences of duplicated busco contigs 

seqtk subseq $pri_asm ${genotype}_duplicated_busco_contigs.txt > ${genotype}_duplicated_busco_contigs.fasta 

## run purge_dup

## Step 1: Run minimap2 to align pacbiodata and generate paf files 

## For PacBio CSS reads 

minimap2 -xasm20 ${genotype}_duplicated_busco_contigs.fasta ${pb_list} -t 26 | gzip -c - > ${genotype}_hifi_reads.busco_dup_contigs.paf.gz
~/apps/purge_dups/bin/pbcstat ${genotype}_hifi_reads.busco_dup_contigs.paf.gz #(produces PB.base.cov and PB.stat files)
~/apps/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log

## Step2: Split genome and do self alignment 
# ~/apps/purge_dups/bin/split_fa ${genotype}_duplicated_busco_contigs.fasta > ${genotype}_duplicated_busco_contigs.split
# minimap2 -xasm5 -DP ${genotype}_duplicated_busco_contigs.split ${genotype}_duplicated_busco_contigs.split -t 28 | gzip -c - > ${genotype}.split.self.busco_dup_contigs.paf.gz

## Step 3: Purge haplotigs and overlaps 
~/apps/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov $${genotype}.split.self.busco_dup_contigs.paf.gz > dups.bed 2> purge_dups.log

## plot coverage 
~/apps/purge_dups/scripts/hist_plot.py PB.base.cov PB.base.png
# ~/apps/purge_dups/scripts/hist_plot.py -c cutoffs PB.base.cov PB.mazia.png

## Step 4: Get purged primary and haplotig sequences from draft assembly
~/apps/purge_dups/bin/get_seqs -e dups.bed ${genotype}_duplicated_busco_contigs.fasta

##~/apps/purge_dups/bin/get_seqs dups.bed $pri_asm > purged.without_e.fa 2> hap.without_e.fa 






















