# reads and parse blastn output (out6 format) of mitochodorial sequences against assembled genomes

library(tidyverse)
library(phylotools)

parse_blast_outformat6 <- function(path_dir="files_directory",
                                   # file_prefix="files_prefix",
                                   file_ending_pattern="files ending pattern including the last separator",
                                   strings_to_remove="strings to remove after the prefix"){
  # lists of blastn files
  
  
  Names = list.files(path = path_dir, pattern=file_ending_pattern) %>% 
    str_remove(strings_to_remove)
  
  # blastp
  
  blast_out = list()
  
  
  for (i in Names) {
    blastp <- read.delim(paste0(path_dir,"/", i,file_ending_pattern), header = F) %>%
      
      dplyr::rename (
        qseqid = V1,
        sseqid = V2,
        qlen = V3,
        slen = V4,
        pident = V5,
        length = V6,
        mismatch = V7,
        gapopen = V8,
        qstart = V9,
        qend = V10,
        sstart = V11,
        send = V12,
        evalue = V13,
        bitscore = V14
      ) %>%
      mutate(quer_sub = i ,
             query_cov = round(length/qlen,1),
             sub_cov = round(length/slen,1)) %>%
      # separate (quer_sub, into = c("query_genome", "subject_genome"), sep="\\." )
      
      blast_out[[i]] <- blastp # add it to your list
      
  }
  
  
  # Unlist  
  
  blast_out = 
    do.call(rbind, blast_out) 
  
  
  # remove row.names()
  
  row.names(blast_out) <- NULL
  #rm(blastp)
  
  # save output as r object 
  
  return(blast_out)
}


blast_out <- parse_blast_outformat6(path_dir = "./", 
                                    file_ending_pattern = "_reads.e99.blastn.out",
                                    strings_to_remove = "refmito.wildb")

blast_out %>%
  filter(pident>80) %>% 
  filter(query_cov>80) %>% 
  select(qseqid) %>% # when the query is Mitochondrial DNA sequences and, query_cov to sub_cov when the reference is EV genome.
  distinct() %>%
  write.table(paste0("refmito_ID", genotype, ".80_80.txt"),
              col.names = F,
              row.names = F,
              quote = F,
              sep = '\t')
