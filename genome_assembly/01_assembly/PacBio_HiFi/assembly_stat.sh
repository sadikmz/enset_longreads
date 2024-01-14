#!/bin/bash

# Convert Hifiasm output in gfa format to fasta, and generate statistics of assembly contiguity and completeness. 
mkdir assembled_fasta summary_assembly_stat

# collect hifiasm output file prefix

ls -lht *bp.*.gfa | awk '{print $9}' | grep -v utg | sed 's/.bp.*//g' | sort | uniq > para_list.txt


param=$(cat para_list.txt)

for i in ${param}

do 

# extract fasta 
    awk '/^S/{print ">"$2"\n"$3}' ${i}.bp.p_ctg.gfa > assembled_fasta/${i}.fasta 
    awk '/^S/{print ">"$2"\n"$3}' ${i}.bp.hap1.p_ctg.gfa > assembled_fasta/${i}.hap1.fasta 
    awk '/^S/{print ">"$2"\n"$3}' ${i}.bp.hap2.p_ctg.gfa > assembled_fasta/${i}.hap2.fasta 

# assembly contiguity 

        ~/apps/gfastats/build/bin/gfastats \
        -f assembled_fasta/${i}.fasta  \
        -t \
        --stats \
        --cmd > summary_assembly_stat/assembly_${i}.stat

        ~/apps/gfastats/build/bin/gfastats \
        -f assembled_fasta/${i}.hap1.fasta  \
        -t \
        --stats \
        --cmd > summary_assembly_stat/assembly_${i}.hap1.stat

        ~/apps/gfastats/build/bin/gfastats \
        -f assembled_fasta/${i}.hap2.fasta  \
        -t \
        --stats \
        --cmd > summary_assembly_stat/assembly_${i}.hap2.stat

# assembly completeness 
busco  -i assembled_fasta/${i}.fasta  -m genome  -l embryophyta_odb10 -c 48 -o busco.${i}.out --long -f --offline 
busco  -i assembled_fasta/${i}.hap1.fasta -m genome -l embryophyta_odb10 -c 48 -o busco.${i}.hap1.out --long -f --offline 
busco  -i assembled_fasta/${i}.hap2.fasta -m genome -l embryophyta_odb10 -c 48 -o busco.${i}.hap2.out --long -f --offline 

done
