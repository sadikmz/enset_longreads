#!/bin/bash
#SBATCH --job-name=mztldk
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


module purge
# module load GCC/8.3.0
# MY_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OMP_NUM_THREADS=$MY_NUM_THREADS

# Telomeres
genotype=ref

for in epo epo_canu mazia wildb wildc 
do 

	# search tidk 
	tidk search EV_${i}_v3.unmasked.fna -s TTTAGGG -o ${i}.T3AG3.tidk --dir tidk_out --extension tsv 
	tidk search EV_${i}_v3.unmasked.fna  -s CCCTAAA -o ${i}.C3TA3.tidk --dir tidk_out --extension tsv 

	## Option #2
	samtools faidx EV_${i}_v3.unmasked.fna  

	cat EV_${i}_v3.unmasked.fna.fai | awk '{print $1, $2}' OFS='\t' > ${i}.genome_size.txt
	bedtools makewindows -g ${i}.genome_size.txt -w 10000 > ${i}_windows.bed 

	bedtools nuc -fi EV_${i}_v3.unmasked.fna   -bed ${i}_windows.bed -pattern TTTAGGG -seq | grep -v '#' | cut -f1,2,3,14 > ${i}.10k.T3AG3.bed
	bedtools nuc -fi EV_${i}_v3.unmasked.fna   -bed ${i}_windows.bed -pattern CCCTAAA -seq | grep -v '#' | cut -f1,2,3,14 > ${i}.10k.C3TA3.bed


	trf EV_${i}_v3.unmasked.fna 2 7 7 80 10 50 500 -f -d -m -h

	python trf2gff -d EV_${i}_v3.unmasked.*.dat -o EV_${i}_trf.gff3

	# Filter redundant information
	cat EV_${i}_trf.gff3 | awk '{split($9,a,";");print $1"\t"$4"\t"$5"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[9]}' > EV_${i}_trf.split.txt


done 


# ## Option #2
# samtools faidx ../${genotype}.fa

# cat ../${genotype}.fa.fai | awk '{print $1, $2}' OFS='\t' > ${genotype}.genome_size.txt
# bedtools makewindows -g ${genotype}.genome_size.txt -w 10000 > ${genotype}_windows.bed 

# bedtools nuc -fi ../${genotype}.fa -bed ${genotype}_windows.bed -pattern TTTAGGG -seq | grep -v '#' | cut -f1,2,3,14 > ${genotype}.10k.T3AG3.bed
# bedtools nuc -fi ../${genotype}.fa -bed ${genotype}_windows.bed -pattern CCCTAAA -seq | grep -v '#' | cut -f1,2,3,14 > ${genotype}.10k.C3TA3.bed


# Centromeres 


# # Install
# git clone https://github.com/oushujun/EDTA.git
# cd EDTA
# conda env create -f EDTA.yml
# conda activate EDTA
# perl EDTA.pl

# # RUN
# source /your/path/anaconda3/bin/activate EDTA

# perl EDTA.pl --genome genome.fa --sensitive 1 --overwrite 1 --anno 1 --species others --threads 10

# grep 'Copia' genome.mod.EDTA.TEanno.gff3 > TE_Copia.split.bed
# grep 'Gypsy' genome.mod.EDTA.TEanno.gff3 > TE_Gypsy.split.bed
# grep 'Helitron' genome.mod.EDTA.TEanno.gff3 > TE_Helitron.split.bed
# grep 'MULE-MuDR' genome.mod.EDTA.TEanno.gff3 > TE_MULE-MuDR.split.bed

## TRF

# Install
conda create -n TRF
conda activate TRF
conda install -c bioconda trf

# RUN (defualt)
# trf yoursequence.fa 2 7 7 80 10 50 500 -f -d -m

# python TRF2GFF.py -d trf_output.dat -o genome_trf.gff3

# # Filter redundant information
# cat genome_trf.gff3 | awk '{split($9,a,";");print $1"\t"$4"\t"$5"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[9]}' > genome_trf.split.txt




grep -w 'period=1' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_1bp.split.bed
grep -w 'period=2' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_2bp.split.bed
grep -w 'period=3' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_3bp.split.bed
grep -w 'period=4' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_4bp.split.bed
grep -w 'period=6' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_6bp.split.bed
grep -w 'period=7' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_7bp.split.bed
grep -w 'period=8' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_8bp.split.bed
grep -w 'period=9' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_9bp.split.bed
grep -w 'period=10' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_10bp.split.bed
grep -w 'period=17' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_17bp.split.bed
grep -w 'period=20' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_20bp.split.bed
grep -w 'period=21' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_21bp.split.bed
grep -w 'period=25' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_25bp.split.bed
grep -w 'period=29' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_29bp.split.bed
grep -w 'period=38' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_38bp.split.bed
grep -w 'period=14' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_14bp.split.bed
grep -w 'period=15' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_15p.split.bed
grep -w 'period=19' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_19bp.split.bed
grep -w 'period=24' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_24bp.split.bed
grep -w 'period=134' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_134bp.split.bed
grep -w 'period=145' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_145bp.split.bed
grep -w 'period=146' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_146bp.split.bed
grep -w 'period=148' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_148bp.split.bed
grep -w 'period=279' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_279bp.split.bed
grep -w 'period=280' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_280bp.split.bed
grep -w 'period=291' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_291bp.split.bed
grep -w 'period=424' EV_mazia_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_mazia_trf_424bp.split.bed
# EV_wildb



grep -w 'period=1' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_1bp.split.bed
grep -w 'period=2' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_2bp.split.bed
grep -w 'period=3' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_3bp.split.bed
grep -w 'period=4' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_4bp.split.bed
grep -w 'period=6' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_6bp.split.bed
grep -w 'period=7' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_7bp.split.bed
grep -w 'period=8' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_8bp.split.bed
grep -w 'period=9' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_9bp.split.bed
grep -w 'period=10' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_10bp.split.bed
grep -w 'period=13' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_13bp.split.bed
grep -w 'period=17' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_17bp.split.bed
grep -w 'period=20' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_20bp.split.bed
grep -w 'period=21' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_21bp.split.bed
grep -w 'period=25' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_25bp.split.bed
grep -w 'period=24' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_24bp.split.bed
grep -w 'period=29' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_29bp.split.bed
grep -w 'period=30' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_30bp.split.bed
grep -w 'period=38' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_38bp.split.bed
grep -w 'period=14' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_14bp.split.bed
grep -w 'period=15' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_15bp.split.bed
grep -w 'period=16' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_16bp.split.bed
grep -w 'period=19' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_19bp.split.bed
grep -w 'period=18' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_18bp.split.bed
grep -w 'period=24' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_24bp.split.bed
grep -w 'period=28' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_28bp.split.bed
grep -w 'period=114' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_114bp.split.bed

grep -w 'period=135' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_135bp.split.bed
grep -w 'period=134' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_134bp.split.bed
grep -w 'period=133' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_133bp.split.bed
grep -w 'period=152' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_152bp.split.bed

grep -w 'period=145' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_145bp.split.bed
grep -w 'period=146' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_146bp.split.bed
grep -w 'period=148' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_148bp.split.bed
grep -w 'period=278' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_278bp.split.bed
grep -w 'period=279' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_279bp.split.bed
grep -w 'period=280' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_280bp.split.bed
grep -w 'period=291' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_291bp.split.bed
grep -w 'period=292' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_292bp.split.bed

grep -w 'period=436' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_436bp.split.bed
grep -w 'period=437' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_437bp.split.bed
grep -w 'period=438' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_438bp.split.bed
grep -w 'period=436' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_436bp.split.bed

grep -w 'period=426' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_426bp.split.bed
grep -w 'period=425' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_425bp.split.bed
grep -w 'period=424' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_424bp.split.bed
grep -w 'period=423' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_423bp.split.bed
grep -w 'period=412' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_412bp.split.bed
grep -w 'period=413' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_413bp.split.bed
grep -w 'period=414' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_414bp.split.bed
grep -w 'period=351' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_351bp.split.bed
grep -w 'period=352' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_352bp.split.bed
grep -w 'period=353' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_353bp.split.bed
grep -w 'period=354' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_354bp.split.bed
grep -w 'period=355' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_355bp.split.bed


grep -w 'period=281' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_281bp.split.bed
grep -w 'period=292' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_292bp.split.bed
grep -w 'period=439' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_439bp.split.bed
grep -w 'period=438' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_438bp.split.bed
grep -w 'period=437' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_437bp.split.bed
grep -w 'period=436' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_436bp.split.bed
grep -w 'period=435' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_435bp.split.bed
grep -w 'period=434' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_434bp.split.bed
grep -w 'period=433' EV_wildb_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildb_trf_433bp.split.bed


# wildc
grep -w 'period=1' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_1bp.split.bed
grep -w 'period=2' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_2bp.split.bed
grep -w 'period=3' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_3bp.split.bed
grep -w 'period=4' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_4bp.split.bed
grep -w 'period=6' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_6bp.split.bed
grep -w 'period=7' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_7bp.split.bed
grep -w 'period=8' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_8bp.split.bed
grep -w 'period=9' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_9bp.split.bed
grep -w 'period=10' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_10bp.split.bed
grep -w 'period=12' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_12bp.split.bed
grep -w 'period=13' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_13bp.split.bed
grep -w 'period=14' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_14bp.split.bed
grep -w 'period=16' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_16bp.split.bed
grep -w 'period=17' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_17bp.split.bed
grep -w 'period=18' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_18bp.split.bed
grep -w 'period=19' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_19bp.split.bed
grep -w 'period=15' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_15bp.split.bed
grep -w 'period=20' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_20bp.split.bed
grep -w 'period=21' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_21bp.split.bed
grep -w 'period=25' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_25bp.split.bed
grep -w 'period=24' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_24bp.split.bed
grep -w 'period=29' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_29bp.split.bed
grep -w 'period=30' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_30bp.split.bed
grep -w 'period=38' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_38bp.split.bed
grep -w 'period=14' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_14bp.split.bed
grep -w 'period=15' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_15bp.split.bed
grep -w 'period=16' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_16bp.split.bed
grep -w 'period=19' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_19bp.split.bed
grep -w 'period=18' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_18bp.split.bed
grep -w 'period=22' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_22bp.split.bed
grep -w 'period=27' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_27bp.split.bed
grep -w 'period=26' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_26bp.split.bed
grep -w 'period=28' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_28bp.split.bed
grep -w 'period=23' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_23bp.split.bed
grep -w 'period=31' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_31bp.split.bed
grep -w 'period=45' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_45bp.split.bed
grep -w 'period=48' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_48bp.split.bed
grep -w 'period=114' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_114bp.split.bed

# grep -w 'period=135' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_135bp.split.bed
grep -w 'period=133' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_133bp.split.bed
grep -w 'period=134' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_134bp.split.bed
grep -w 'period=133' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_133bp.split.bed
grep -w 'period=152' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_152bp.split.bed
grep -w 'period=135' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_135bp.split.bed
grep -w 'period=426' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_426bp.split.bed
grep -w 'period=157' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_157bp.split.bed
grep -w 'period=148' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_148bp.split.bed


grep -w 'period=145' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_145bp.split.bed
grep -w 'period=146' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_146bp.split.bed
grep -w 'period=144' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_144bp.split.bed
grep -w 'period=152' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_152bp.split.bed
grep -w 'period=278' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_278bp.split.bed
grep -w 'period=279' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_279bp.split.bed
grep -w 'period=280' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_280bp.split.bed
grep -w 'period=291' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_291bp.split.bed
grep -w 'period=292' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_292bp.split.bed

grep -w 'period=351' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_351bp.split.bed

grep -w 'period=436' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_436bp.split.bed
grep -w 'period=437' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_437bp.split.bed
grep -w 'period=438' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_438bp.split.bed
grep -w 'period=436' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_436bp.split.bed

grep -w 'period=426' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_426bp.split.bed
grep -w 'period=425' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_425bp.split.bed
grep -w 'period=424' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_424bp.split.bed
grep -w 'period=423' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_423bp.split.bed
grep -w 'period=412' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_412bp.split.bed
grep -w 'period=413' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_413bp.split.bed
grep -w 'period=414' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_414bp.split.bed
grep -w 'period=351' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_351bp.split.bed
grep -w 'period=352' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_352bp.split.bed
grep -w 'period=353' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_353bp.split.bed
grep -w 'period=354' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_354bp.split.bed
grep -w 'period=355' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_355bp.split.bed


grep -w 'period=281' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_281bp.split.bed
grep -w 'period=292' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_292bp.split.bed
grep -w 'period=413' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_413bp.split.bed
grep -w 'period=424' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_424bp.split.bed

grep -w 'period=439' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_439bp.split.bed
grep -w 'period=438' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_438bp.split.bed
grep -w 'period=437' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_437bp.split.bed
grep -w 'period=436' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_436bp.split.bed
grep -w 'period=435' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_435bp.split.bed
grep -w 'period=434' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_434bp.split.bed
grep -w 'period=433' EV_wildc_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_wildc_trf_433bp.split.bed


## Epo

grep -w 'period=1' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_1bp.split.bed
grep -w 'period=2' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_2bp.split.bed
grep -w 'period=3' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_3bp.split.bed
grep -w 'period=4' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_4bp.split.bed
grep -w 'period=6' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_6bp.split.bed
grep -w 'period=7' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_7bp.split.bed
grep -w 'period=8' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_8bp.split.bed
grep -w 'period=9' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_9bp.split.bed
grep -w 'period=10' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_10bp.split.bed
grep -w 'period=12' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_12bp.split.bed
grep -w 'period=13' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_13bp.split.bed
grep -w 'period=14' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_14bp.split.bed
grep -w 'period=16' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_16bp.split.bed
grep -w 'period=17' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_17bp.split.bed
grep -w 'period=18' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_18bp.split.bed
grep -w 'period=19' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_19bp.split.bed
grep -w 'period=15' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_15bp.split.bed
grep -w 'period=20' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_20bp.split.bed
grep -w 'period=21' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_21bp.split.bed
grep -w 'period=25' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_25bp.split.bed
grep -w 'period=24' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_24bp.split.bed
grep -w 'period=29' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_29bp.split.bed
grep -w 'period=30' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_30bp.split.bed
grep -w 'period=38' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_38bp.split.bed
grep -w 'period=14' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_14bp.split.bed
grep -w 'period=15' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_15bp.split.bed
grep -w 'period=16' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_16bp.split.bed
grep -w 'period=19' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_19bp.split.bed
grep -w 'period=18' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_18bp.split.bed
grep -w 'period=22' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_22bp.split.bed
grep -w 'period=27' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_27bp.split.bed
grep -w 'period=26' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_26bp.split.bed
grep -w 'period=28' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_28bp.split.bed
grep -w 'period=23' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_23bp.split.bed
grep -w 'period=31' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_31bp.split.bed
grep -w 'period=45' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_45bp.split.bed
grep -w 'period=48' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_48bp.split.bed
grep -w 'period=114' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_114bp.split.bed

# grep -w 'period=135' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_135bp.split.bed
grep -w 'period=133' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_133bp.split.bed
grep -w 'period=134' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_134bp.split.bed
grep -w 'period=133' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_133bp.split.bed
grep -w 'period=152' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_152bp.split.bed
grep -w 'period=135' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_135bp.split.bed
grep -w 'period=426' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_426bp.split.bed
grep -w 'period=157' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_157bp.split.bed
grep -w 'period=148' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_148bp.split.bed


grep -w 'period=145' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_145bp.split.bed
grep -w 'period=146' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_146bp.split.bed
grep -w 'period=144' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_144bp.split.bed
grep -w 'period=152' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_152bp.split.bed
grep -w 'period=278' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_278bp.split.bed
grep -w 'period=279' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_279bp.split.bed
grep -w 'period=280' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_280bp.split.bed
grep -w 'period=291' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_291bp.split.bed
grep -w 'period=292' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_292bp.split.bed

grep -w 'period=351' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_351bp.split.bed

grep -w 'period=436' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_436bp.split.bed
grep -w 'period=437' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_437bp.split.bed
grep -w 'period=438' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_438bp.split.bed
grep -w 'period=436' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_436bp.split.bed

grep -w 'period=426' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_426bp.split.bed
grep -w 'period=425' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_425bp.split.bed
grep -w 'period=424' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_424bp.split.bed
grep -w 'period=423' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_423bp.split.bed
grep -w 'period=412' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_412bp.split.bed
grep -w 'period=413' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_413bp.split.bed
grep -w 'period=414' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_414bp.split.bed
grep -w 'period=351' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_351bp.split.bed
grep -w 'period=352' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_352bp.split.bed
grep -w 'period=353' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_353bp.split.bed
grep -w 'period=354' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_354bp.split.bed
grep -w 'period=355' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_355bp.split.bed


grep -w 'period=281' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_281bp.split.bed
grep -w 'period=292' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_292bp.split.bed
grep -w 'period=413' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_413bp.split.bed
grep -w 'period=424' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_424bp.split.bed

grep -w 'period=439' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_439bp.split.bed
grep -w 'period=438' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_438bp.split.bed
grep -w 'period=437' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_437bp.split.bed
grep -w 'period=436' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_436bp.split.bed
grep -w 'period=435' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_435bp.split.bed
grep -w 'period=434' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_434bp.split.bed
grep -w 'period=433' EV_epo_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_trf_433bp.split.bed
# epo_canu
grep -w 'period=1' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_1bp.split.bed
grep -w 'period=2' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_2bp.split.bed
grep -w 'period=3' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_3bp.split.bed
grep -w 'period=4' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_4bp.split.bed
grep -w 'period=6' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_6bp.split.bed
grep -w 'period=7' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_7bp.split.bed
grep -w 'period=8' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_8bp.split.bed
grep -w 'period=9' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_9bp.split.bed
grep -w 'period=10' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_10bp.split.bed
grep -w 'period=12' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_12bp.split.bed
grep -w 'period=13' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_13bp.split.bed
grep -w 'period=14' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_14bp.split.bed
grep -w 'period=16' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_16bp.split.bed
grep -w 'period=17' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_17bp.split.bed
grep -w 'period=18' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_18bp.split.bed
grep -w 'period=19' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_19bp.split.bed
grep -w 'period=15' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_15bp.split.bed
grep -w 'period=20' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_20bp.split.bed
grep -w 'period=21' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_21bp.split.bed
grep -w 'period=25' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_25bp.split.bed
grep -w 'period=24' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_24bp.split.bed
grep -w 'period=29' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_29bp.split.bed
grep -w 'period=30' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_30bp.split.bed
grep -w 'period=38' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_38bp.split.bed
grep -w 'period=14' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_14bp.split.bed
grep -w 'period=15' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_15bp.split.bed
grep -w 'period=16' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_16bp.split.bed
grep -w 'period=19' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_19bp.split.bed
grep -w 'period=18' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_18bp.split.bed
grep -w 'period=22' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_22bp.split.bed
grep -w 'period=27' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_27bp.split.bed
grep -w 'period=26' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_26bp.split.bed
grep -w 'period=28' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_28bp.split.bed
grep -w 'period=23' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_23bp.split.bed
grep -w 'period=31' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_31bp.split.bed
grep -w 'period=45' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_45bp.split.bed
grep -w 'period=48' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_48bp.split.bed
grep -w 'period=114' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_114bp.split.bed

# grep -w 'period=135' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_135bp.split.bed
grep -w 'period=133' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_133bp.split.bed
grep -w 'period=134' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_134bp.split.bed
grep -w 'period=133' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_133bp.split.bed
grep -w 'period=152' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_152bp.split.bed
grep -w 'period=135' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_135bp.split.bed
grep -w 'period=426' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_426bp.split.bed
grep -w 'period=157' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_157bp.split.bed
grep -w 'period=148' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_148bp.split.bed


grep -w 'period=145' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_145bp.split.bed
grep -w 'period=146' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_146bp.split.bed
grep -w 'period=144' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_144bp.split.bed
grep -w 'period=152' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_152bp.split.bed
grep -w 'period=278' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_278bp.split.bed
grep -w 'period=279' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_279bp.split.bed
grep -w 'period=280' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_280bp.split.bed
grep -w 'period=291' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_291bp.split.bed
grep -w 'period=292' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_292bp.split.bed

grep -w 'period=351' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_351bp.split.bed

grep -w 'period=436' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_436bp.split.bed
grep -w 'period=437' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_437bp.split.bed
grep -w 'period=438' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_438bp.split.bed
grep -w 'period=436' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_436bp.split.bed

grep -w 'period=426' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_426bp.split.bed
grep -w 'period=425' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_425bp.split.bed
grep -w 'period=424' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_424bp.split.bed
grep -w 'period=423' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_423bp.split.bed
grep -w 'period=412' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_412bp.split.bed
grep -w 'period=413' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_413bp.split.bed
grep -w 'period=414' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_414bp.split.bed
grep -w 'period=351' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_351bp.split.bed
grep -w 'period=352' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_352bp.split.bed
grep -w 'period=353' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_353bp.split.bed
grep -w 'period=354' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_354bp.split.bed
grep -w 'period=355' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_355bp.split.bed


grep -w 'period=281' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_281bp.split.bed
grep -w 'period=292' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_292bp.split.bed
grep -w 'period=413' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_413bp.split.bed
grep -w 'period=424' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_424bp.split.bed

grep -w 'period=439' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_439bp.split.bed
grep -w 'period=438' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_438bp.split.bed
grep -w 'period=437' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_437bp.split.bed
grep -w 'period=436' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_436bp.split.bed
grep -w 'period=435' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_435bp.split.bed
grep -w 'period=434' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_434bp.split.bed
grep -w 'period=433' EV_epo_canu_trf.split.txt | sed 's/copies=//g' | cut -f1-3,6 > EV_epo_canu_trf_433bp.split.bed


