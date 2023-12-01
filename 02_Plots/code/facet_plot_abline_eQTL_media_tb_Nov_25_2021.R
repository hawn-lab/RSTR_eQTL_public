#Call Library
library(tidyverse)
source("publication/plot_eQTL.R")

#set run time
run_time <- gsub(":","_", Sys.Date())

#Data
load("data_clean/eQTL_data_all.RData")
# snp annotation
snp.annot <- read_tsv("data/snp/snp.annot_megaEX_unfiltered_n_80_filtered_by_annot_201118.txt")

#Call eQTL Genes
eQTL_final <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv") %>% 
  select(gene, snps, rsID) %>% 
  rename(snp=snps)

#Make plots
plot_ls <- plot_eQTL(to_plot=eQTL_final, snp_anno=snp.annot, 
                     type="cis", rename_rsID=TRUE)

# Save
for(i in 1:length(plot_ls)){
  print(i)
  filename <- paste0("02_Plots/00_Media_TB_Plots/", run_time, "_", 
                     names(plot_ls)[i],".pdf")
  ggsave(filename=filename, plot_ls[[i]], width=4, height=4)
}
