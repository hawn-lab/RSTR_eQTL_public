library(tidyverse)
library(patchwork)
source("publication/plot_eQTL.R")
source("publication/plot_eQTL_delta.R")
#
#### Data ####
#Call eQTL Genes
eQTL_final <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv") %>% 
  rename(snp=snps)
eQTL_lead <- eQTL_final %>% filter(lead_SNP=="Y")

load("data_clean/eQTL_data_all.RData")

# snp annotation
snp.annot <- read_tsv("data/snp/snp.annot_megaEX_unfiltered_n_80_filtered_by_annot_201118.txt")

# LM results
cis_eqtl <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv") %>% 
  mutate(status="Media") %>% 
  bind_rows(read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv") %>% 
              mutate(status="+Mtb")) %>% 
  rename(snp=snps) %>% 
  select(snp, gene, status,FDR) %>% 
  pivot_wider(names_from = status, values_from = FDR) %>% 
  mutate(across(c(Media, `+Mtb`), 
                ~ifelse(is.na(.), "FDR > 0.1",
                        paste("FDR =", signif(., digits = 2)))))

trans_eqtl <- read_csv("01_Cytokine/results/delta_model/2023-06-23_cytokine_interaction_co.csv") %>% 
  rename(snp=snps) %>% 
  select(snp, gene, p.value)

#### Figure 3 ####
#Define SNP and gene for plotting
genes3 <- c("SIRPB1","CPSF1","GPX2","GSTO2","MLYCD","PRKAG2")
to_plot3 <- eQTL_lead %>% 
  filter(gene %in% genes3) %>% 
  select(gene, snp, rsID) %>% 
  #Add FDR
  left_join(cis_eqtl) %>% 
  mutate(gene = factor(gene, levels=genes3))  %>% 
  arrange(gene, snp)

plot_ls3 <- plot_eQTL(to_plot=to_plot3, snp_anno=snp.annot, 
                      type="cis", rename_rsID = TRUE, FDR = TRUE)

plot3 <- plot_spacer() + plot_spacer() + plot_ls3[[1]] + plot_ls3[[2]] + 
  plot_ls3[[3]] + plot_ls3[[4]] + plot_ls3[[5]] + plot_ls3[[6]] +
  plot_annotation(tag_levels = list(c("b","c","d","e","f","g")), tag_suffix = ".") +
  plot_layout(guides = 'collect', ncol=4)
plot3

ggsave(file="publication/archive/Fig3_eQTL.express.pdf", plot3, 
       width = 12.2, height = 9)

#### Figure 4 ####
#Define SNP and gene for plotting
to_plot4 <- read_csv("01_Cytokine/results/delta_model/2023-06-23_cytokine_interaction_co.csv") %>% 
  filter(p.value < 0.05) %>% 
  select(gene, snps, rsID) %>% 
  rename(snp=snps) %>% 
  filter(snp !="rs55966428") %>% 
  mutate(gene = factor(gene, levels=c("IFNB1","TNF","IL1B","IL6"))) %>% 
  arrange(gene, snp) %>% 
  left_join(trans_eqtl)

plot_ls4 <- plot_eQTL_delta(to_plot=to_plot4, snp_anno=snp.annot, cis_anno=eQTL_final,
                      type="trans", rename_rsID = TRUE, FDR = TRUE)

plot4 <- wrap_plots(plot_ls4) +  
  plot_annotation(tag_levels = "a", tag_suffix = ".") +
  plot_layout(guides = 'collect', ncol=4)
plot4

ggsave(filename = "publication/Fig4_eQTL.cyto.pdf", plot4, width=12, height=9)
ggsave(filename = "publication/Fig4_eQTL.cyto.png", plot4, width=12, height=9)

#### Figure 5? ####
#BMP6
to_plot5a <- eQTL_lead %>% 
  filter(gene %in% c("BMP6")) %>% 
  select(gene, snp, rsID) %>% 
  #Add FDR
  left_join(cis_eqtl) %>% 
  arrange(gene, snp)
  
to_plot5b <- read_csv("01_Cytokine/results/delta_model/2023-06-23_cytokine_interaction_co.csv") %>% 
  filter(p.value < 0.05) %>% 
  select(gene, snps, rsID, p.value) %>% 
  rename(snp=snps) %>% 
  filter(snp =="rs55966428") %>% 
  arrange(gene, snp)

plot_ls5a <- plot_eQTL(to_plot=to_plot5a, snp_anno=snp.annot, cis_anno=eQTL_final,
                      type="cis", rename_rsID = TRUE, FDR=TRUE)

plot_ls5b <- plot_eQTL_delta(to_plot=to_plot5b, snp_anno=snp.annot, cis_anno=eQTL_final,
                       type="trans", rename_rsID = TRUE, FDR=TRUE)

plot5 <- plot_ls5a$BMP6_rs55966428 + plot_ls5b$IFNB1_rs55966428 +  
  plot_annotation(tag_levels = "a", tag_suffix = ".") +
  plot_layout(guides = 'collect', widths = c(2,1))
# plot5

ggsave(filename = "publication/Fig5_BMP6.pdf", plot5, width=6.1, height=4.5)
ggsave(filename = "publication/Fig5_BMP6.png", plot5, width=6.1, height=4.5)

#### Figure S1 ####
# all eQTL
to_plots1 <- eQTL_lead %>% 
  select(gene, snp, rsID) %>% 
  filter(!gene %in% c(genes3,"BMP6")) %>% 
  #Add FDR
  left_join(cis_eqtl) %>% 
  arrange(gene, snp)

plot_lss1 <- plot_eQTL(to_plot=to_plots1, snp_anno=snp.annot, cis_anno=eQTL_final,
                      type="cis", rename_rsID = TRUE, FDR=TRUE)

plots1 <- wrap_plots(plot_lss1) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") +
  plot_layout(guides = 'collect', ncol=5)
# plots1

ggsave(plots1, filename = "publication/FigS1.all.eQTL.pdf", 
       width=15, height=9)
# ggsave(plots1, filename = "publication/FigS1.all.eQTL.png", 
#        width=15, height=9)


#### Figure SX. media and tb cyto (instead of delta) ####
to_plotsx <- read_csv("01_Cytokine/results/delta_model/2023-06-23_cytokine_interaction_co.csv") %>% 
  filter(p.value < 0.05) %>% 
  select(gene, snps, rsID) %>% 
  rename(snp=snps) %>% 
  mutate(gene = factor(gene, levels=c("IFNB1","TNF","IL1B","IL6"))) %>% 
  arrange(gene, snp) %>% 
  left_join(trans_eqtl)

plot_lssx <- plot_eQTL(to_plot=to_plotsx, snp_anno=snp.annot, cis_anno=eQTL_final,
                            type="trans", rename_rsID = TRUE, FDR = TRUE)

plotSX <- wrap_plots(plot_lssx) +  
  plot_annotation(tag_levels = "a", tag_suffix = ".") +
  plot_layout(guides = 'collect', ncol=3)
plotSX

ggsave(plotSX, filename = "publication/FigSX_eQTL.cyto.notDelta.pdf", 
       width=10, height=12)
