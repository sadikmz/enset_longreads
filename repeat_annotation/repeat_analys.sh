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
path_dir=absolute_path_to_pseudochromosome_assemblies
genotype= #mazia/wildb/wildc/epo

## Run RepeatProteinMask
RepeatProteinMask ${genotype}.final.fna.00.s1k.0.fa -noLowSimple -pvalue 0.0001 -engine ncbi

## RepeatMasker 

RepeatMasker -species Viridiplantae ${genotype}.fna -nolow -no_is -norna -engine ncbi -pa 12 -gff -poly -small -excln 
mkdir rmbed_out
RM2Bed.py --out_dir rmbed_out/ ${genotype}.fna.out 
bedtools getfasta -fi ${genotype}.fna -bed rmbed_out/${genotype}.fna_rm.bed -fo ${genotype}.RM.fa 


## Tandem Repeat Finde, TRF, Piler, RepeatModeler2
trf ${genotype}.fna 2 5 7 80 10 50 2000 -m -h 
trf2gff -o - < ${genotype}.fna.*.dat > ${genotype}.trf.gff3

## RepeatModeler 
BuildDatabase -name wildc ${genotype}
RepeatModeler -database wildc -pa 27 -LTRStruct -genotypeSampleSizeMax 480MB -threads 48 

# # Miniature Inverted-repeat Transposable Elements (MITE)-Hunter 
# ~/apps/MITE-hunter/MITE_Hunter_manager.pl -i ${genotype}.final.fna -g ${genotype}_MH -c 28 -n 5 -S 12345678 –P 1
# ## MITE_Hunter_manager TEs:  ${genotype}_MH_Step8_singlet.fa 

## combine all 

cat ${genotype}.final.fna.RM.fa ${genotype}.piler_library.fa ${genotype}.families.fa ${genotype}.trf.fa ${genotype}.RPM.fa > ${genotype}.redun_repeat_lib.fa 

# Remove redundant sequences
cd-hit -i ${genotype}.redun_repeat_lib.fa  -d 0 -o ${genotype}.non_redun_repeat_lib.fa  -c 0.80 -n 5 -G 1 -g 0 -M 0 -T 24

#
rep_file=${genotype}.non_redun_repeat_lib.fa 

#makeblastdb -in /home/lifesci/lfrwtp/data/musa.ensete.ensemblmocot.db.v1/musa.ensete.ensemblmocot.fa -dbtype prot

blastx \
      -db /home/lifesci/lfrwtp/data/musa.ensete.ensemblmocot.db.v1/musa.ensete.ensemblmocot.fa \
      -query $rep_file \
      -outfmt '6 qseqid qlen slen qstart qend sstart send stitle pident length evalue bitscore qcovhsp scovhsp' \
      -evalue 1e-10 \
      -max_hsps 1 \
      -num_threads 24 \
      -out ${rep_file}.blastx.out

##### convert blast hits to bed file and exclude the bed file interval from fasta

cut -f 1,4,5,8 ${rep_file}.blastx.out |\
# remove any transposable element hits
grep -v -E 'polyprotein|Polyprotein|Retrotransposon|retro|Retro|RNA-directed|transposable|Transposable|gypsy|Gypsy|Copia|copia|DNA/RNA polymerases|Integrase' |\
awk '{if ($2>$3) print $1,$3,$2,".",".","-"; else print $1,$2,$3,".",".","+";}' OFS='\t' | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ${rep_file}.protein_cleared.draft.bed
bedtools merge -i ${rep_file}.protein_cleared.draft.bed -s -c 6 -o distinct > ${rep_file}.protein_cleared.final.bed
# generate the index file
samtools faidx $rep_file
# reformate the index file to genome file
cut -f1,2 ${rep_file}.fai > ${rep_file}.txt
# final bed file to exclude
bedtools complement -i ${rep_file}.protein_cleared.final.bed -g $rep_file.txt > ${rep_file}.final.bed
# final protein coding gene free final de novo repeats
bedtools getfasta -fi $rep_file -bed ${rep_file}.final.bed > $rep_file.final

cut -f 1,4,5 ${rep_file}.blastx.out | awk '{if ($2>$3) print $1,$3,$2,".",".","-"; else print $1,$2,$3,".",".","+";}' OFS='\t' | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ${rep_file}.protein_cleared.draft.bed
bedtools getfasta -fi -bed ${rep_file}.protein_cleared.draft.bed > ${genotype}.non_redun_repeat_lib_prots_excluded.fa 

## filter out sequence < 80 bases 

./filter_fasta_length.py ${rep_file}.final 80 > ${genotype}.non_redun_repeat_lib_prots_excluded.80.fa 

## Run EDTA

genome=${genotype}.fna
type=all
species=others
threads=28
rmout=path_to_repeatmasker_${genotype}.fna.out
cds=EV_mazia_bedadeti.cds.fa # CDS from shorted assembeld enset genomes 

EDTA.pl \
-genome $genome \
--species $species \
--overwrite 0 \
--sensitive 1 \
--anno 1 \
--evaluate 1 \
--threads $threads \
--rmout $rmout \
--cds $cds \
--step anno

## Run RepeatMasker 

## repeat library 
cat ${genotype}.non_redun_repeat_lib_prots_excluded.80.fa EDTA_out/${genotype}.fna.mod.EDTA.TElib.fa > ${genotype}.repeat_lib.fa

lib=${genotype}.repeat_lib.fa

RepeatMasker  -pa 12 -lib $lib -cutoff 225 -gff ${genotype}.fna -poly -a -u -xsmall

# Estimates of repeatitive seqeunce divergence 
calcDivergenceFromAlign.pl -s ${genotype}.divsum  ${genotype}.fna.align
createRepeatLandscape.pl -div ${genotype}.divsum > ${genotype}.repeatlandscape.html
createRepeatLandscape.pl -div ${genotype}.divsum  -t "Mazia" > ${genotype}.title_added.repeatlandscape.html


## extract unknown or unspecified repeats
rmsk2bed < ${genotype}.fna.out | cut -f1-3,11 | grep 'Unspecified\|Unknown' | cut -f1-3 | sort -V | uniq > ${genotype}unknown_unspecified.RM.out.bed








