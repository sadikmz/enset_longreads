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
  - [NLR literature and benchmarking](#nlr-literature-and-benchmarking)
    - [Assembly quality check](#assembly-quality-check)
    - [Population structure](#population-structure)

## Organnels genome assembly (Chloroplast and mitochondorial genomes)
- Chloroplast genome: excluding reads before assembly 
After self-correction and low-quality region trimming, SMRT long reads were mapped to a previously released patchouli chloroplast genome (NCBI accession: NC_042796.1)56 by Minimap2 (version 2.5-r572)57 with the parameters -t 96 -ax map-pb. As the mapped reads belonged to chloroplasts, the remaining unmapped reads were first assembled to obtain the nuclear genome <https://www.nature.com/articles/s41467-022-31121-w#Sec10>
- <https://doi.org/10.1111/1755-0998.13616>



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
HiC -> HiC-Pro -> Juicerbox -> 3D-DNA -> HiC-plotter (ONT) <https://www.nature.com/articles/s42003-022-03646-9#Sec11>
HiFi -> HiC -> HiC -> 3D-DNA -> Juicebox manual inspection <https://www.frontiersin.org/articles/10.3389/fpls.2022.1012277/full#h3>
CLR -> HiC -> LACHESIS <https://www.mdpi.com/2079-7737/11/10/1492>

### Generating and vizualing HiC contact probability maps

- 3D-DNA and Juicebox <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08697-0#Sec12>
- [HiCPlotter](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0767-1)

### Benchmarking genome articles
Chinese pepper  <https://doi.org/10.1111/pbi.13926>
Allium crops <https://www.nature.com/articles/s41467-022-34491-3>
Tea genome <https://doi.org/10.1038/s41438-020-0288-2>
Turmeric <https://www.frontiersin.org/articles/10.3389/fpls.2022.1003835/full>
Sesame <https://doi.org/10.1016/j.jare.2022.10.004>
Nibea coibor /HiFi+HiC/ <https://doi.org/10.1038/s41597-022-01804-6>
Arabidopsis /Nanopore + HiFi/ <https://doi.org/10.1016/j.gpb.2021.08.003>
Velvet bean <https://doi.org/10.1093/dnares/dsac031>
<https://onlinelibrary.wiley.com/doi/epdf/10.1111/tpj.15968>
Watermelon <https://doi.org/10.1093/molbev/msac168>
Brassicaceae diploid Orychophragmus violaceus <https://doi.org/10.1016/j.xplc.2022.100432>
Poales species Carex cristatella <https://doi.org/10.1093/g3journal/jkac211>

## NLR literature and benchmarking
<https://doi.org/10.1016/j.pbi.2022.102311>
Insight into the structure and molecular mode of action of plant paired NLR immune receptors  <https://doi.org/10.1042/EBC20210079>
Receptor-mediated nonhost resistance in plants  <https://doi.org/10.1042/EBC20210080>
Unconventional R proteins in the botanical tribe Triticeae  <https://doi.org/10.1042/EBC20210081>
Evolution of resistance (R) gene specificity <https://doi.org/10.1042/EBC20210077>
Activation and Regulation of NLR Immune receptors <https://doi.org/10.1093/pcp/pcac116>
Tsw – A case study on structure-function puzzles in plant NLRs with unusually large LRR domains <https://www.frontiersin.org/articles/10.3389/fpls.2022.983693/full>
### Assembly quality check 

### Population structure
<https://www.frontiersin.org/articles/10.3389/fpls.2022.1022169/full#h3>
