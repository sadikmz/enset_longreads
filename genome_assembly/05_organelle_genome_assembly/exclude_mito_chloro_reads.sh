#!/bin/bash
#SBATCH --job-name=wbmito using picard
#SBATCH --partition=cnodeICARD}/picard.jar SortSam \
#SBATCH --nodes=12_out/"$base"_refmito.sam   \
#SBATCH --ntasks-per-node=1se"_refmito.sorted.sam   \
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



# specify the directory containing the input HiFi reads files and the reference genome
PICARD=/home/lifesci/lfrwtp/miniconda3/envs/picard/share/picard-2.27.5-0
reads_dir="/home/lifesci/lfrwtp/data/hifi_data/"
reference="mito_chlo.db/mitoRefSeq_MA_PD_EV_chrpt.fa"
genotype= #mazia/wildb/wild/epo

# loop through all fastq files in the directory
for file in "$reads_dir"/*.adapt.discarded.fasta; do
## step1: map HiFi reads against mitochorial sequences of NCBI mito_refeseq MA and PA?
    # remove the path and extension to get the base name
    base=$(basename "$file" | sed 's/.adapt.discarded.fasta//g')
    echo "Processing $base"

    # map the HiFi reads to the reference genome using minimap2
    minimap2 -xasm20 "$reference" "$file" -t 28 -a -o minimap2_out/"$base"_refmito.sam

    # # sort the SAM file using picard
    java -Xmx50G -jar ${PICARD}/picard.jar SortSam \
        I=minimap2_out/"$base"_refmito.sam   \
        O=minimap2_out/"$base"_refmito.sorted.sam   \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        TMP_DIR=./tmp/
        VALIDATION_STRINGENCY=SILENT


        # extract the mapped reads using samtools

    # # convert the BAM file to SAM file using picard
    java -Xmx50G -jar $PICARD/picard.jar SamFormatConverter \
        I=minimap2_out/"$base"_refmito.sorted.mapped.sam \
        O=minimap2_out/"$base"_refmito.sorted.mapped.bam

    # step3: Extract sequeces on HiFi reads mapped with mito sequeces, and excluded the mapped mito sequences from HiFi reads  
    # generated bed file 
    bedtools bamtobed -i minimap2_out/"$base"_refmito.sorted.mapped.bam | bedtools sort | bedtools merge | awk '{print $1}' | sort | uniq > minimap2_out/"$file"_mapped_reads.mitochrpt.txt 
    bedtools bamtobed -i minimap2_out/"$base"_refmito.sorted.mapped.bam | bedtools sort | bedtools merge > minimap2_out/"$file"_mapped_reads.mitochrpt.bed  

    # get reads name and exclude mito_chrpt mapped reads  
    samtools faidx ${path_dir}/"$file".adapt.discarded.fasta | awk '{print $1}' | sort | uniq > minimap2_out/"$file".reads_name.txt
    grep -v -f minimap2_out/"$file"_mapped_reads.mitochrpt.txt minimap2_out/"$file".reads_name.txt > minimap2_out/"$file"_unmapped_reads.txt   

    # get fasta for unmapped reads 
    seqtk subseq ${path_dir}/"$file".adapt.discarded.fasta minimap2_out/"$file"_unmapped_reads.txt > minimap2_out/"$file".adapt_mitochrpt_discarded.fasta

    # or 
    # remove mito sequence from reads 
    awk '{print $1 $2}' OFS='\t' ${path_dir}/"$file".adapt.discarded.fasta  > ${path_dir}/"$file".genome_size.txt  
    bedtools complement -i minimap2_out/"$file"_refmito.reads.bed -g ${path_dir}/"$file".genome_size.txt > minimap2_out/"$file"_mitochrpt_excluded.bed
    # final protein coding gene free final de novo repeats
    bedtools getfasta -fi ${path_dir}/"$file".adapt.discarded.fasta -bed minimap2_out/"$file".mitochrpt_excluded.bed > minimap2_out/"$file".adapt_mitochrpt_discarded.complement.fasta

done 
