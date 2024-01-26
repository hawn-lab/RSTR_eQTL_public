library(tidyverse)

#### Table S1. cis eQTL ####
# eQTL in MEDIA condition
eQTL_MEDIA <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv",
                       col_types = cols(
                         snps = col_character(),
                         gene = col_character(),
                         statistic = col_double(),
                         pvalue = col_double(),
                         FDR = col_double(),
                         beta = col_double()))%>%
  select(snps, gene, beta, FDR) %>% 
  rename(estimate_media=beta, FDR_media=FDR)

# eQTL in TB condition
eQTL_TB <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv",
                    col_types = cols(
                      snps = col_character(),
                      gene = col_character(),
                      statistic = col_double(),
                      pvalue = col_double(),
                      FDR = col_double(),
                      beta = col_double()))%>%
  select(snps, gene, beta, FDR) %>% 
  rename(estimate_mtb=beta, FDR_mtb=FDR)

# eQTL with significant interaction lm
interact_lm <- read_csv("00_eQTL/results/interaction_model/2023-06-21_eQTL_interaction_co.csv") %>%
  select(snps, gene, Estimate, FDR) %>% 
  rename(estimate_interact=Estimate, FDR_interact=FDR)

eQTL_all <- full_join(eQTL_TB, eQTL_MEDIA) %>% 
  full_join(interact_lm) %>% 
  #add position data
  left_join(read_tsv("data/snp/snpsloc_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt"),
            by=c("snps"="id")) %>% 
  #add rsID
  left_join(read_csv("data/snp/eQTL_rsid_key.csv"),
            by=c("snps"="SNPS")) %>% 
  select(snps, rsID, chr, pos, everything()) %>% 
  mutate(rsID = ifelse(rsID==snps, NA, rsID)) %>% 
  arrange(gene, snps)

write_csv(eQTL_all, file="publication/TableS1.cis.eQTL.csv")

#### Table S5. cytokine eQTL ####
eQTL_cyto <- read_csv("01_Cytokine/results/interaction_model/2023-06-21_cytokine_interaction_co.csv") %>% 
  select(snps,rsID, chr, pos, gene, Estimate, FDR) %>% 
  rename(estimate_interact=Estimate, FDR_interact=FDR) %>% 
  mutate(rsID = ifelse(rsID==snps, NA, rsID)) %>% 
  arrange(gene, snps)

write_csv(eQTL_cyto, file="publication/TableS5.cyto.eQTL.csv")

#### Table S3. enrichment ####
load("00_eQTL/results/enrichment.RData")
map <- read_csv("publication/rrvgo_map.csv") %>% 
  select(term, parentTerm) %>% 
  rename(parentPathway=parentTerm)

enrich_all <- bind_rows(enrich_c2, enrich_c5) %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>% 
  mutate(term = gsub("GOBP_","",pathway),
         term = gsub("_"," ", tolower(term))) %>% 
  left_join(map) %>% 
  select(gs_cat, gs_subcat, pathway, parentPathway, group_in_pathway, 
         size_pathway,`k/K`, FDR, genes) %>% 
  mutate(genes = as.character(genes))

write_csv(enrich_all, file="publication/TableS3.eQTL.enrichment.csv")
