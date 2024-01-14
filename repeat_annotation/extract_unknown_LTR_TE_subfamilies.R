## add TEsorter annotation of EDTA_intact_LTR_TE from reformat_tesorter_intact_LTR_TE


file_path="/Users/u1866313/rstudio/pacbio_clr_css/longreads_project/repeat_annotation/EDTA_ltr_age_estimatation_updated/"

remove_strings=".intact_ltr_TE_tesorter_annot.intersect.txt"
file_pattern=".intact_ltr_TE_tesorter_annot.intersect.txt"


intact_ltr_TE_tesorter_annot <- list()


for (i in list.files(path = file_path,pattern = file_pattern)){
  
  # capture basename 
  BASENAME <- basename(i) %>%
    str_remove(remove_strings)
  
  # reads TEsoter output 
  
  intact_TEsorter_annotation <- 
  read.delim(paste0(file_path,i), header = F) %>% 
    mutate(genome = BASENAME) 
    intact_ltr_TE_tesorter_annot[[i]] <- intact_TEsorter_annotation
}


intact_ltr_TE_tesorter_annot = 
  do.call(rbind, intact_ltr_TE_tesorter_annot) 

rownames(intact_ltr_TE_tesorter_annot) <- NULL

intact_ltr_TE_tesorter_annot %>%
  head()

## combine EDTA TEs with their TEsorter annotation 

# EDTA_intact_LTR_TE_tesorter_annot_combined <-
  EDTA_intact_LTR_TE %>% 
  # filter(genome == "epo") %>%
  mutate(TSD_5P_start = as.numeric(TSD_5P_start),
         TSD_3P_end = as.numeric(TSD_3P_end)) %>%
  left_join(
    intact_ltr_TE_tesorter_annot %>%
      dplyr::rename(TSD_5P_start = V2,
                    TSD_3P_end = V3)) %>%
  # filter(str_detect(TSD_5P_start,"10378989")) %>%
  # str()
  left_join(
    mazia_dante_unknown %>%
      separate(seqid,c("seqid","TSD_5P_start","TSD_3P_end")) %>%
      separate(B,c("superfam","ltr_fam_unk"),sep = "\\|") %>%
      mutate(
        TSD_5P_start = as.numeric(TSD_5P_start),
        TSD_3P_end = as.numeric(TSD_3P_end),
        superfam = str_replace(superfam,"/","_")))  %>%
  left_join(
    EV_ltr_age_combined %>%
      rename(LTR_3P_start_RT=LTR_5P_start_RT, 
             LTR_3P_end_RT=LTR_5P_end_RT) %>%
      mutate(LTR_3P_start_RT = as.numeric(LTR_3P_start_RT))
  ) %>%
  # select(genome,superfamily,TE_type,V7,V8,V9,insertion_time,TSD_5P_start,TSD_3P_end, superfam, ltr_fam_unk, ltr_age_mya, gnom) %>%
  mutate(len = TSD_3P_end - TSD_5P_start +1 ,
         insertion_time = round(insertion_time/1000000,2),
         len = round(len/1000,2),
         TE_type = str_replace_na(TE_type,"NA"),
         V7 = str_replace_na(V7,"NA"),
         V7 = str_replace(V7,"\\.",""),
         V8 = str_replace(V8,"\\.",""),
         V9 = str_replace(V9,"\\.",""),
         TE_type = case_when(str_detect(V7,"Ty1") ~ "Ty1_copia",
                             str_detect(V7,"Ty3") ~ "Ty3_gypsy",
                             TRUE ~ TE_type),
         TE_type = case_when(TE_type=="NA" ~ superfamily,
                             TRUE ~ TE_type),
         TE_type = str_replace(TE_type,"unknown","Unknown"),
         ltr_fam= str_c(V8,V9,sep = "/"),
         ltr_fam = case_when(ltr_fam=="/" ~ "Unknown",
                             TRUE ~ ltr_fam),
         superfamily = case_when(str_detect(V7,"gypsy") ~ "Gypsy",
                                 str_detect(V7,"copia") ~ "Copia",
                                 TRUE ~ superfamily),
         superfamily = str_replace(superfamily,"unknown","Unknown"),
         genome = case_when(genome=="epo" ~ "Epo",
                            genome=="mazia" ~ "Mazia",
                            genome=="wildb" ~ "Wild-B",
                            genome=="wildc" ~ "Wild-C"),
         ltr_fam = str_replace_na(ltr_fam,"NA"),
         ltr_fam = case_when(ltr_fam=="Ivana/" ~ "Ivana",
                             ltr_fam=="Tork/" ~ "Tork",
                             ltr_fam=="Ikeros/" ~ "Ikeros",
                             ltr_fam=="Angela/" ~ "Angela",
                             ltr_fam=="Bianca/" ~ "Bianca",
                             ltr_fam=="Ale/" ~ "Ale",
                             ltr_fam=="Bryco/" ~ "Bryco",
                             ltr_fam=="SIRE/" ~ "SIRE",
                             ltr_fam=="Lyco/" ~ "Lyco",
                             ltr_fam=="Gymco-I/" ~ "Lyco",
                             ltr_fam=="Lyco/" ~ "Lyco",
                             ltr_fam=="Gymco-III/" ~ "Gymco-III",
                             ltr_fam=="TAR/" ~ "TAR",
                             ltr_fam=="Osser/" ~ "Osser",
                             ltr_fam=="Gymco-II/" ~ "Gymco-II",
                             ltr_fam=="Ty1-outgroup/" ~ "Ty1-outgroup",
                             # ltr_fam=="NA" && superfamily == "LTR" ~ "LTR/unknown",
                             TRUE ~ ltr_fam 
         ),         
         TE_type = case_when(TE_type=="LTR" &  superfamily=="Gypsy" ~ "Ty3_gypsy",
                             TE_type=="LTR" &  superfamily=="Copia" ~ "Ty1_copia",
                             # ltr_fam=="NA" && superfamily == "LTR" ~ "LTR/unknown",
                             TRUE ~ TE_type
         ),         
  ) %>%
  # select(genome,superfamily,TE_type,ltr_fam,insertion_time,len,superfam, ltr_fam_unk, ltr_age_mya, gnom) %>%
  mutate(
    ltr_fam = case_when(ltr_fam=="NA" ~ ltr_fam_unk,
                        TRUE ~ ltr_fam)
  ) %>%
  dplyr::rename(
    Superfamily = TE_type,
    Superfamily0 = superfamily
  ) %>% 
  filter(len > 0) %>%
  mutate(
    ltr_fam = str_replace(ltr_fam,"non-chromovirus/OTA","non-chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/Galadriel","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/CRM","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/Tcn1","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"non-chromovirus/non-chromo-outgroup","non-chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/Tekay","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/chromo-unclass","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/chromo-outgroup","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/chromo-unclass","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"non-chromovirus/Phygy","non-chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/Chlamyvir","chromovirus"),
    ltr_fam = str_replace(ltr_fam,"chromovirus/Reina","chromovirus"),
    ltr_fam = str_replace_na(ltr_fam, "Unknown"),
    ltr_fam_unk = str_replace_na(ltr_fam_unk,"NA"),
    Superfamily = case_when(Superfamily=="Unknown" & ltr_fam=="Tork" ~ "Ty1_copia",
                            Superfamily=="Unknown" & ltr_fam=="chromovirus" ~ "Ty3_gypsy",
                            Superfamily=="Unknown" & ltr_fam=="non-chromovirus" ~ "Ty3_gypsy",
                            Superfamily=="LTR" & superfam=="Ty1_copia" ~ "Ty1_copia",
                            Superfamily=="LTR" & superfam=="Ty3_gypsy" ~ "Ty3_gypsy",
                            TRUE ~ Superfamily),
    
    
    
  )   %>%
    select(genome,seqid,TSD_5P_start, TSD_3P_end,Superfamily,ltr_fam,insertion_time,len,superfam, ltr_fam_unk, ltr_age_mya, gnom) %>%
    filter(ltr_fam == "Unknown") %>%
    filter(genome == "Mazia") %>%
    select(seqid,TSD_5P_start,TSD_3P_end) %>%
    write.table(paste0(file_path,"mazia_unknonw_LTR_TE.bed"),
                col.names = F, row.names = F, quote = F, sep = '\t')