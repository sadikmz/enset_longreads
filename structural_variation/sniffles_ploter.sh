#!/bin/bash
#SBATCH --job-name=ngmlrmz
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=fc.%J.out
#SBATCH --error=fc.%J.err
#SBATCH --mail-type=FAIL,END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=Sadik.Muzemil@warwick.ac.uk


module purge

python3 sniffles2_plots.py -i ngmlr_sniffle_vcf
# python3 sniffles2_plots.py -i <file_name> -o <output_folder>

# parse sniffles vcf 
python3 sniffles2_plots.py -i ngmlr_sniffle_vcf
cat ngmlr_sniffle_vcf/wildb_ngmlr.picard.vcf | python sniffles2_vcf_parser.py parsesv  > ngmlr_sniffle_vcf/wildb_ngmlr_sniffles.txt 
cat ngmlr_sniffle_vcf/wildc_ngmlr.picard.vcf | python sniffles2_vcf_parser.py parsesv  > ngmlr_sniffle_vcf/wildc_ngmlr_sniffles.txt 
cat ngmlr_sniffle_vcf/mazia_ngmlr.picard.vcf | python sniffles2_vcf_parser.py parsesv  > ngmlr_sniffle_vcf/mazia_ngmlr_sniffles.txt 
cat ngmlr_sniffle_vcf/epo_ngmlr.picard.vcf | python sniffles2_vcf_parser.py parsesv  > ngmlr_sniffle_vcf/epo_ngmlr_sniffles.txt 
