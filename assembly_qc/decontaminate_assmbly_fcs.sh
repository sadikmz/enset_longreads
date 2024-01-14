!/bin/bash

# [FCS-adapter](https://github.com/ncbi/fcs) was used for decontaminated foriegn genomes sequences in enset genomes. 

SHM_LOC=/shared/software/ncbi

data_dir=data_assembly
outdir=decontam
taxid=4639


assembled_genomes=$(ls -lht ${data_dir}.fna | awk '{print $9)')

for i in ${assembled_genomes}; do

	# base name of assembeld genomes
	BASENAME=$(basename $i | sed 's/adapt_plastome_discarded.fna//g')
	
	# run 
	python3 ${SHM_LOC}/run_fcsgx.py \
	--fasta ${BASENAME}.adapt_plastome_discarded.fna \
	--out-dir ${BASENAME}.outputdir \
	--gx-db "${SHM_LOC}/gxdb/all" \
	--tax-id 4639
done 

# Discard contaminant sequences 

## extract coordinates of foreign sequence
cat ${genome}.outputdir/${genome}.adapt_plastome_discarded.4639.fcs_gx_report.txt | grep -v '#' | awk '{print $1, $2, $3}' OFS='\t' > ${genome}.outputdir/${genome}_contaminants.bed 
## contaminats contigs id and whole sequence 
cat ${genome}.outputdir/${genome}.adapt_plastome_discarded.4639.fcs_gx_report.txt | grep -v '#' | awk '{print $1}' > ${genome}.outputdir/${genome}_contaminants_contigs.txt 
seqtk subseq ${genome}.adapt_plastome_discarded.fna ${genome}.outputdir/${genome}_contaminants_contigs.txt  > ${genome}.outputdir/${genome}_contaminants_contigs.fasta 

## break contigs before and after the contaminant regions 
samtools faidx ${genome}_contaminants_contigs.fasta 
cat ${genome}_contaminants_contigs.fasta.fai | awk '{print $1, $2}' > ${genome}.outputdir/${genome}_contaminants_contigs.size
## final bed file to exclude
bedtools complement -i ${genome}.outputdir/${genome}_contaminants.bed -g ${genome}.outputdir/${genome}_contaminants_contigs.size > ${genome}.outputdir/${genome}_contaminants_excluded.bed 
## final contaminats excluded contigs region 
bedtools getfasta -fi ${genome}.outputdir/${genome}_contaminants_contigs.fasta -bed ${genome}.outputdir/${genome}_contaminants_excluded.bed > ${genome}.outputdir/${genome}_contaminants_excluded_contigs.fasta  

## add contaminants excluded fasta into the assembly and remove temporary files 
samtools faidx ${genome}.outputdir/${genome}.adapt_plastome_discarded.fna 
cat ${genome}.adapt_plastome_discarded.fna.fai | grep -v '#' | awk '{print $1}' | grep -v -f ${genome}.outputdir/${genome}_contaminants_contigs.txt  > ${genome}.outputdir/${genome}_decontaminated_contigs.txt
seqtk subseq ${genome}.adapt_plastome_discarded.fna ${genome}.outputdir/${genome}_decontaminated_contigs.txt > ${genome}.outputdir/${genome}_decontaminated_contigs.fasta  
cat ${genome}.outputdir/${genome}_decontaminated_contigs.fasta ${genome}.outputdir/${genome}_contaminants_excluded_contigs.fasta > ${genome}.outputdir/${genome}.adapt_plastome_discarded_decontam.fasta

## manual check contigs lenght to discard shorter contigs  
samtools faidx ${genome}.outputdir/${genome}.adapt_plastome_discarded_decontam.fasta
cat ${genome}.outputdir/${genome}.adapt_plastome_discarded_decontam.fasta.fai | awk '{print $1,$2}' | sort -k2 -n > ${genome}.outputdir/sorted.${genome}_contaminants_excluded_contigs.size

## remove temp files  
rm ${genome}.outputdir/${genome}_decontaminated_contigs.fasta ${genome}.outputdir/${genome}_contaminants_excluded_contigs.fasta 
