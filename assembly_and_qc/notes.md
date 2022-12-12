# Merqury plot for PacBio HiFi assemblies 

## spectral plots


<img src="merqury.mazia_hifi.out.mazia_s33_adapt_discarded.spectra-cn.st.png" title="Mazia specra-cn plot" width="40%"> <img src="merqury.wild_c_s33_adapt_discarded.spectra-cn.st.png" title="Wild-C specra-cn plot" width="40%"> 

<img src="merqury.wild_b.out.wild_b_s33_adapt_discarded.spectra-cn.st.png" title="Wild-B specra-cn plot" width="40%">


The spectral plot shows that genome contains very diverse haplotypes see <https://github.com/marbl/merqury/issues/59>

Modifying purge_dup cutoff appears to make no significant difference base don the discussion here <https://github.com/dfguan/purge_dups/issues/14>

<img src="merqury.mazia_hifi.out.purged.spectra-cn.st.png" title="Mazia after modifying purge_dup cutoff" width="40%">

cutoff modified: 2 7 11 13 20 39
cutoff initial: 5 7 11 13 22 39


hamming error rate - measures teh occurrence of the k-mers from the unexpected haplotypes to the other in the entire assembly see <https://github.com/marbl/merqury/issues/53> \n

Also there are edge cases in HiFi assemblies where HiFi or Illumina sequencing biases are involved (e.g. homopolymer / 2-mer microsatellite indel errors, GC biases), so having 100% completeness seems very suspicious... <https://github.com/marbl/merqury/issues/84>

Hi-C scaffolding <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05087-x>
and <https://github.com/zhaotao1987/Plant_genome_Hi-C_scaffolding>

Genome <https://plantmethods.biomedcentral.com/articles/10.1186/s13007-022-00964-1>

Missing wild ancestro - banana <https://www.frontiersin.org/articles/10.3389/fpls.2022.969220/full>

Annotation liftoff and other <https://www.biorxiv.org/content/10.1101/2022.12.01.518658v1.full>
and <https://www.biorxiv.org/content/10.1101/2022.12.01.518724v1.full>

Dotplot <https://github.com/WHops/nahrchainer>
Sequence identify hitmap centromer etc <https://eichlerlab.gs.washington.edu/help/glogsdon/Shared_with_Pille/StainedGlass_adjustedScale.R>
