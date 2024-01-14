# PILER as described in https://www.drive5.com/piler

genotype= #mazia/wildb/wildc/epo
genome=${genotype}.fna 

## Step 1.  Local alignment:

### 1.1 Split the genome into chunks small enough for PALS
pyfasta split -n 10 ${genome}.final.fna

### Collect splitted chunks as fofn
for i in ${genome}.final.fna.??
  do echo ${i} >> file_list.txt
  done

### Test pals with one the chunks

### 1.2  align each chunk to itself using the -self option of PALS
files=$(cat file_list.txt)

for i in ${files}
   do
     pals -self ${i} -out ${i}.gff
   done

###  1.3 Align each different pair of chunks to each other using the ‑query and ‑target options of PALS
for i in ${files}
   do
     pals -query ${i} -target ${i} -out ${i}_tq.hits.gff
   done

### 1.3 Concatenate the hit files produced in steps 1.1 and 1.2 into a single hit file.
# cat ../mazia_genomic_v1.??.fas.masked.masked_tq.hits.gff ../mazia_genomic_v1.??.fas.masked.masked.hits.gff > piler_self_tq.hits.gff
cat *.gff > piler_self_tq.hits.gff

## Step 2. TRS (Transposed Repeat Signature) search
# use output gff file from step as input here

piler -trs piler_self_tq.hits.gff -out ${genome}_trs.gff -piles coordinate_${genome}_trs.gff

## Step 3: Library construction
#Library construction takes the genome sequence and TRS output file as input, and produces a library of consensus sequences as output.
#3a. Making a family FASTA file

mkdir fams
piler -trs2fasta ${genome}_trs.gff -seq ${genome}.final.fna -path fams
mkdir aligned_fams

#3b. Making a family alignment
cd fams
for fam in *
    do
        muscle -in $fam -out ../aligned_fams/$fam -maxiters 1 -diags1
    done

cd ..
mkdir cons
cd aligned_fams
for fam in *
    do
       piler -cons $fam -out ../cons/$fam -label $fam
    done

cd ../cons
cat * > ../piler_library.fasta
