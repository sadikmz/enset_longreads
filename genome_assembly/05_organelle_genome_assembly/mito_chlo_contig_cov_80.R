library(tidyverse)
genotype= #mazia/wildb/wildc/epo
read.delim(paste0("blast_out/", genotype, ".blastn.sorted.merged.bed"), header=F) %>% 
left_join(read.delim(paste0("blast_out/", genotype, ".contig.sizes"), header=F) %>% 
                rename (V4=V2 ))%>% 
         mutate ( cov = (V3-V2)/V4) %>% filter (cov > 0.80 )  %>% 
select(V1) %>% 
distinct() %>% write.table(paste0("blast_out/", genotype, "_mito_chrpt_contigs.cov_80.txt"), 
	col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')