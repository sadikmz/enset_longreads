library(ape)
mazia_dante_unknown_unspecific <-
  read.gff("~/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/dante_unknown_unspecific/epo_uknown_unspecific_dante/domains_filtered.gff") %>%
  mutate(genome = "epo") %>%
  rbind(
    read.gff("~/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/dante_unknown_unspecific/mazia_uknown_unspecific_dante/domains_filtered.gff") %>%
      mutate(genome = "mazia"),
    read.gff("~/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/dante_unknown_unspecific/wildb_uknown_unspecific_dante/domains_filtered.gff") %>%
      mutate(genome = "wildb"),
    read.gff("~/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/dante_unknown_unspecific/wildc_uknown_unspecific_dante/domains_filtered.gff") %>%
      mutate(genome = "wildc")
  ) %>%
  separate(attributes,into = c("A","B","C","D","E","F","G","H","I","J","K","L","M"), sep = ";" )%>%
  select(seqid,genome,B) %>%
  mutate(B = str_remove(B,"Final_Classification=Class_\\w+\\|LTR\\|")) %>%
  mutate(B = str_remove_all(B,"\\|OTA\\|Tat\\|Retand")) %>%
  mutate(B = str_remove_all(B,"\\|OTA\\|Tat")) %>%
  mutate(B = str_remove_all(B,"\\|Reina")) %>%
  mutate(B = str_remove_all(B,"\\|Galadriel")) %>%
  mutate(B = str_remove(B,"Final_Classification\\=Class_II\\|")) %>%
  mutate(B = str_remove(B,"Final_Classification\\=Class_I|")) %>%
    # mutate(B = str_replace(B,"chromovirus\\/Reina","chromovirus")) %>%
  # mutate(B = str_replace(B,"chromovirus\\/CRM","chromovirus")) %>%
  # mutate(B = str_replace(B,"chromovirus\\/Tcn1","chromovirus")) %>%
  # mutate(B = str_replace(B,"non-chromovirus\\/non-chromo-outgrou","non-chromovirus")) %>%
  # mutate(B = str_replace(B,"chromovirus\\/Tekay","chromovirus")) %>%
  # mutate(B = str_replace(B,"chromovirus\\/Chlamyvir","chromovirus")) %>%
  # mutate(B = str_replace(B,"chromovirus\\/chromo-outgroup","chromovirus")) %>%
  # mutate(B = str_replace(B,"non-chromovirus/OTA","non-chromovirus")) %>%
  separate(seqid,into = c("seqid","coord"), sep = ":") %>%
  separate(coord, into = c("start","end"), sep = "-") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(end),
           len = end - start+1) %>%
    distinct()
    
mazia_dante_unknown_unspecific %>%
#   filter(len == "NA")
# filter(genome == "wildc") %>%
#     nrow()
  group_by(genome,B) %>%
    summarise(count = n(), tot_len = sum(len)) %>%
    as.data.frame() %>%
    arrange(desc(count))
  distinct(B) 


mazia_dante_unknown %>%
  distinct(B)

read.gff("~/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/EDTA_ltr_age_estimatation_updated/dante_unknown/wildc_uknown_dante/domains_filtered.gff") %>%
  separate(attributes,into = c("A","B","C","D","E","F","G","H","I","J","K","L","M"), sep = ";" )%>%
  select(seqid,B) %>%
  mutate(B = str_remove(B,"Final_Classification=Class_\\w+\\|LTR\\|")) %>%
  mutate(B = str_remove_all(B,"\\|OTA\\|Tat\\|Retand")) %>%
  mutate(B = str_remove_all(B,"\\|OTA\\|Tat")) %>%
  mutate(B = str_remove_all(B,"\\|Reina")) %>%
  mutate(B = str_remove_all(B,"\\|Galadriel")) %>%
  distinct() %>%
  head()
