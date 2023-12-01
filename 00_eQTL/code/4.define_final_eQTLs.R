library(tidyverse)
#Date 
run_time <- gsub(":","_", Sys.Date())

# eQTL in MEDIA condition
eQTL_MEDIA <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv",
                    col_types = cols(
                      snps = col_character(),
                      gene = col_character(),
                      statistic = col_double(),
                      pvalue = col_double(),
                      FDR = col_double(),
                      beta = col_double()))%>%
  dplyr::filter(FDR < 0.01) %>%
  select(snps, gene, FDR) %>% 
  rename(FDR_media=FDR)

eQTL_MEDIA_keep <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv",
                       col_types = cols(
                         snps = col_character(),
                         gene = col_character(),
                         statistic = col_double(),
                         pvalue = col_double(),
                         FDR = col_double(),
                         beta = col_double()))%>%
  select(snps, gene, FDR) %>% 
  rename(FDR_media=FDR)

# eQTL in TB condition
eQTL_TB <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv",
                 col_types = cols(
                   snps = col_character(),
                   gene = col_character(),
                   statistic = col_double(),
                   pvalue = col_double(),
                   FDR = col_double(),
                   beta = col_double()))%>%
  dplyr::filter(FDR < 0.01) %>%
  select(snps, gene, FDR) %>% 
  rename(FDR_mtb=FDR)

# eQTL with significant interaction lm
interact_lm <- read_csv("00_eQTL/results/interaction_model/2023-06-21_eQTL_interaction_co.csv") %>%
  filter(FDR < 0.1) %>% 
  select(snps, gene, FDR) %>% 
  rename(FDR_interact=FDR)

# eQTL are defined as significant in TB condition as well as for interaction term but non-signif in media condition
eQTL_final <- inner_join(eQTL_TB, interact_lm) %>% 
  anti_join(eQTL_MEDIA) %>% 
  #Add media FDR
  left_join(eQTL_MEDIA_keep) %>% 
  #add position data
  left_join(read_tsv("data/snp/snpsloc_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt"),
            by=c("snps"="id")) %>% 
  #add rsID
  left_join(read_csv("data/snp/eQTL_rsid_key.csv"),
            by="snps")

#Tag leads SNP
eQTL_final <- eQTL_final %>% 
  group_by(gene) %>% 
  slice_min(FDR_mtb, n=1, with_ties = TRUE) %>% 
  ungroup() %>% 
  select(snps, gene) %>% 
  mutate(lead_SNP="Y") %>% 
  full_join(eQTL_final)

write_csv(eQTL_final, 
          file=paste0("00_eQTL/results/", run_time, "_final_eQTL_definition.csv"))
