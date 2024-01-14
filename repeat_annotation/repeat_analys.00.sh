#!/bin/bash
#SBATCH --job-name=RM_mz
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --gres=gpu:quadro_rtx_6000:3
#SBATCH --mem-per-cpu=3700
#SBATCH --time=12:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


module purge
module load GCC/8.3.0
#MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$MY_NUM_THREAD

#module load DMTCP/2.6.0

# path to final assembly 
path_dir=/home/lifesci/lfrwtp/hifi_assembly/post_rapid_curation/final_asm/


#module load DMTCP/2.6.0
#dmtcp_launch -i 3600 ./script


#genotype=EV_epo
## Run RepeatProteinMask
#RepeatProteinMask ${genotype}.final.fna.00.s1k.0.fa -noLowSimple -pvalue 0.0001 -engine ncbi

## combine all 

genome=EV_mazia

#cat ${genome}.final.fna.RM.fa ${genome}.piler_library.fa ${genome}.families.fa ${genome}.trf.fa ${genome}.RPM.fa > ${genome}.redun_repeat_lib.fa 

# Remove redundant sequences
#cd-hit -i ${genome}.redun_repeat_lib.fa  -d 0 -o ${genome}.non_redun_repeat_lib.fa  -c 0.80 -n 5 -G 1 -g 0 -M 0 -T 24

#
rep_file=${genome}.non_redun_repeat_lib.fa 

#makeblastdb -in /home/lifesci/lfrwtp/data/musa.ensete.ensemblmocot.db.v1/musa.ensete.ensemblmocot.fa -dbtype prot

#blastx \
#       -db /home/lifesci/lfrwtp/data/musa.ensete.ensemblmocot.db.v1/musa.ensete.ensemblmocot.fa \
#       -query $rep_file \
#       -outfmt '6 qseqid qlen slen qstart qend sstart send stitle pident length evalue bitscore qcovhsp scovhsp' \
#       -evalue 1e-10 \
#       -max_hsps 1 \
#       -num_threads 24 \
#       -out $rep_file.blastx.out

##### convert blast hits to bed file and exclude the bed file interval from fasta
#for i in ${files}
#   do
#      cut -f 1,4,5,8 $rep_file.blastx.out |\
      ## remove any transposable element hits
#      grep -v -E 'polyprotein|Polyprotein|Retrotransposon|retro|Retro|RNA-directed|transposable|Transposable|gypsy|Gypsy|Copia|copia|DNA/RNA polymerases|Integrase' |\
#      awk '{if ($2>$3) print $1,$3,$2,".",".","-"; else print $1,$2,$3,".",".","+";}' OFS='\t' | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > $rep_file.protein_cleared.draft.bed
#      bedtools merge -i $rep_file.protein_cleared.draft.bed -s -c 6 -o distinct > $rep_file.protein_cleared.final.bed
      # generate the index file
#      samtools faidx $rep_file
      # reformate the index file to genome file
#      cut -f1,2 $rep_file.fai > $rep_file.txt
      # final bed file to exclude
#      bedtools complement -i $rep_file.protein_cleared.final.bed -g $rep_file.txt > $rep_file.final.bed
      # final protein coding gene free final de novo repeats
#      bedtools getfasta -fi $rep_file -bed $rep_file.final.bed > $rep_file.final
#   done
####

#cut -f 1,4,5 ${rep_file}.blastx.out | awk '{if ($2>$3) print $1,$3,$2,".",".","-"; else print $1,$2,$3,".",".","+";}' OFS='\t' | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > $rep_file.protein_cleared.draft.bed
#bedtools getfasta -fi -bed $rep_file.protein_cleared.draft.bed > ${genome}.non_redun_repeat_lib_prots_excluded.fa 


## filter out sequence < 80 bases 

#./filter_fasta_length.py $rep_file.final 80 > ${genome}.non_redun_repeat_lib_prots_excluded.80.fa 

## Run RepeatMasker 
genome=EV_mazia

## repeat library 
cat ../repeatmasker_lib/${genome}.non_redun_repeat_lib_prots_excluded.80.fa ../EDTA/${genome}.final.v1.fna.mod.EDTA.TElib.fa > ${genome}.repeat_lib.fa

lib=${genome}.repeat_lib.fa

RepeatMasker -pa 48 -lib $lib -cutoff 225 -gff $genome.final.v1.fna -poly