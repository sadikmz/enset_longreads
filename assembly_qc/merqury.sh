#!/bin/bash
#SBATCH --job-name=merq
#SBATCH --partition=cnode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=00:30:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk

# files/path_to_files 
genoype= #mazia/wildb/wildc/epo
MERQURY=/home/lifesci/lfrwtp/apps/merqury
PATH_WGS=/home/lifesci/lfrwtp/data/ev_wgs
pri_asm=data/${genotype}_s33_k61_adapt_discarded_hic.fasta
pacbio_hifi=/home/lifesci/lfrwtp/working_dir/hifi_assembly_analaysis/organelles_genome/mitochodoria/mito_chlro_free_hifi_reads/mazia_hifi.mito_chlr.excluded.adapt.discarded.fasta.gz
## kmer 
k=20
threads=28

## Get the right k size

## generate meryl kmer database 

$MERQURY/best_k.sh 59000000 > ${genome}.best_kmer.txt 

# create meryl db
meryl k=$k count ${PATH_WGS}/*fq.gz $pacbio_hifi output ${genome}.meryl threads=$threads

# combine meryl db
meryl union-sum output ${genome}.db.meryl ${genome}.meryl threads=$threads

# run merqury for quiver_pilon gap filled assembly 
$MERQURY/merqury.sh ${genome}.db.meryl ${pri_asm} merqury.${genome}.out

# heterozygocity
meryl histogram ${genome}.db.meryl > ${genome}.heterozygocity.hist

# Meryl-based QV

$MERQURY/eval/qv.sh ${genome}.db.meryl ${pacbio_hifi} meryl_qv
