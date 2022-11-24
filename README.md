# PacBio CLR and HiFi reads assembly, annotation, and analysis 

<details open="open">
<summary>Table of Contents</summary>

- [PacBio CLR and HiFi reads assembly, annotation, and analysis](#pacbio-clr-and-hifi-reads-assembly-annotation-and-analysis)
  - [Organnels genome assembly (Chloroplast and mitochondorial genomes)](#organnels-genome-assembly-chloroplast-and-mitochondorial-genomes)
  - [PacBio CLR assembly and scaffolding](#pacbio-clr-assembly-and-scaffolding)
  - [PacBio HiFi assembly and scaffolding](#pacbio-hifi-assembly-and-scaffolding)
    - [Assembly](#assembly)
    - [scaffolding](#scaffolding)
    - [Generating and vizualing HiC contact probability maps](#generating-and-vizualing-hic-contact-probability-maps)
    - [Benchmarking genome articles](#benchmarking-genome-articles)
    - [Assembly quality check](#assembly-quality-check)

## Organnels genome assembly (Chloroplast and mitochondorial genomes)
- Chloroplast genome: excluding reads before assembly 
After self-correction and low-quality region trimming, SMRT long reads were mapped to a previously released patchouli chloroplast genome (NCBI accession: NC_042796.1)56 by Minimap2 (version 2.5-r572)57 with the parameters -t 96 -ax map-pb. As the mapped reads belonged to chloroplasts, the remaining unmapped reads were first assembled to obtain the nuclear genome <https://www.nature.com/articles/s41467-022-31121-w#Sec10>

## PacBio CLR assembly and scaffolding

Quality check on HiFi reads

GeneScope.FK <https://github.com/thegenemyers/GENESCOPE.FK>

PloidyPlot <https://github.com/thegenemyers/MERQURY.FK>

Asset <https://github.com/dfguan/asset#pg>
source <https://academic.oup.com/dnaresearch/article/29/6/dsac035/6815630#379962482>

## PacBio HiFi assembly and scaffolding

### Assembly
Hifiasm with default followed by purge-dups <https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03707-5#Sec1>


### scaffolding 
- Align with BWA/Bowtie (-end-to-end, -very-sensitive -L 30)<https://www.nature.com/articles/s41467-022-34206-8#Sec9>
- HiC-Pro <https://www.nature.com/articles/s41597-022-01671-1#Sec2>
- Hi-CUP <https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-022-01388-y#Sec10>
- LACHESIS
- ALLHiC
- Juicebox after LACHESIS and filtering HiC reads (PacBioCLR) <https://www.nature.com/articles/s42003-022-04145-7#Sec10>

HiC -> LACHESIS

### Generating and vizualing HiC contact probability maps

- 3D-DNA and Juicebox <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08697-0#Sec12>
- [HiCPlotter](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0767-1)

### Benchmarking genome articles

Allium crops <https://www.nature.com/articles/s41467-022-34491-3>

### Assembly quality check 