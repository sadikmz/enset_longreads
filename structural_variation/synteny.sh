#!/bin/bash
#SBATCH --job-name=satsuma
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=3850
#SBATCH --time=48:00:00
#SBATCH --account=su007-mg
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


module purge
module load GCC/10.2.0 OpenMPI/4.0.5

ln -s /home/l/lfrwtp/repeat_annotation/unmasked_reoriented/chrs_only/EV_mazia.fna mazia.fna
ln -s /home/l/lfrwtp/repeat_annotation/unmasked_reoriented/chrs_only/EV_wildb.fna wildc.fna

ref_genome=mazia.fna
query_genome=wildb.fna
outdir=mazRef_wildbqry
mkdir mazRef_wildbqry
threads=128


SatsumaSynteny \
-q $ref_genome \
-t $query_genome \
-n $threads \
-dups true \
-ni 124 \
-n 123 \
-m 1 \
-o $outdir 

