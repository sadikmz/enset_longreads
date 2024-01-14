mkdir jellyfish_histo
dir=/home/l/lfrwtp/data/hifi_cleaned_reads


mkdir jellyfish_histo
dir=/home/l/lfrwtp/data/hifi_cleaned_reads
genotype= #mazia/wildb/wildc/ 

for k in 17 21 27 31 37 41 47 51 57 61; 

do 
	# Generate Kmer histogram
	~/apps/jellyfish-linux count -C -m $k -s 560M --bf-size 50G -t 128 <(zcat ${dir}/${genotype}.mito_chlr.excluded.adapt.discarded.fasta.gz) -o wildc_${k}_mer.reads.jf
	~/apps/jellyfish-linux histo -t 24 wildc_${k}_mer.reads.jf > jellyfish_histo/wildc_${k}_mer.histo
	rm -rf wildc_${k}_mer.reads.jf

	# Estimates of genome hetetrozygousity 
	genomescope2 \
	-i jellyfish_histo/wildc_${k}_mer.histo \
	-o wildc_${i}_gs \
	-n wildc_${k}_mer
	-p 2 \
	-k ${k} \
	--testing \
	--fitted_hist \
	--start_shift \
	--typical_error
done
