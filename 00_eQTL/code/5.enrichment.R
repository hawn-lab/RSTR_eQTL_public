library(tidyverse)
library(SEARchways)

eQTL_final <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv")
genes.OI <- unique(eQTL_final$gene)

# enrich_h <- BIGprofiler(list("eQTL"=genes.OI), category = "H")
enrich_c2 <- BIGprofiler(list("eQTL"=genes.OI), category = "C2", subcategory = "CP")
enrich_c5 <- BIGprofiler(list("eQTL"=genes.OI), category = "C5", subcategory = "GO:BP")


enrich_all <- bind_rows(enrich_c2 %>% filter(FDR < 0.1),
                        enrich_c5 %>% filter(FDR < 0.01))

save(enrich_c2, enrich_c5, file="00_eQTL/results/enrichment.RData")
