library(tidyverse)
library(BIGpicture)
# library(patchwork)
library(rrvgo)
library(GO.db)

#### Data ####
eQTL_final <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv")
genes.OI <- unique(eQTL_final$gene)

#### Format enrichment ####
attach("00_eQTL/results/enrichment.RData")
enrich_c2_signif <- enrich_c2 %>% filter(FDR < 0.1) %>% 
  mutate(pathway = gsub("_"," ",pathway),
         pathway = tolower(pathway),
         pathway = gsub("reactome","REACTOME",pathway),
         pathway = gsub("kegg","KEGG",pathway))

# Group GO terms
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`) %>% 
  mutate(GO_path = tolower(GO_path),
         GO_path = gsub("-"," ", GO_path))

c5_format <- enrich_c5 %>% 
  filter(FDR < 0.1 & group_in_pathway > 1) %>%
  mutate(GO_path = gsub("GOBP_","",pathway),
         GO_path = gsub("_"," ",GO_path),
         GO_path = tolower(GO_path)) %>%
  left_join(goterms)

# Group like terms based on gene content.
simMatrix <- calculateSimMatrix(c5_format$GOID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.9,
                                orgdb="org.Hs.eg.db")
unique(reducedTerms$parentTerm) %>% length()
# scatterPlot(simMatrix, reducedTerms)

##Extract data
x <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)

df <- cbind(as.data.frame(x$points),
            reducedTerms[match(rownames(x$points), reducedTerms$go),
                         c("term", "parent", "parentTerm", "size")]) 
write_csv(df, "publication/rrvgo_map.csv")

df2 <- df %>% 
  mutate(term = gsub("-"," ", term),
         term = tolower(term)) %>% 
  left_join(c5_format, by=c("term"="GO_path")) %>% 
  mutate(term = paste("GOBP", term)) %>% 
  unnest(genes) %>% 
  group_by(parentTerm) %>% 
  dplyr::summarise(genes = list(unique(genes))) %>% 
  ungroup() %>% 
  mutate(pathway = paste("GOBP",parentTerm),
         FDR = 0,
         group_in_pathway = Inf,
         `k/K`=1)

# "GOBP fatty acid metabolic process"="GOBP fatty acid metabolic process\nGOBP regulation of lipid metabolic process\nGOBP positive regulation of peptidyl-tyrosine phosphorylation",
# "GOBP mitotic metaphase plate congression"="GOBP mitotic metaphase plate congression\nGOBP mitotic spindle organization"
enrich_to_plot <- bind_rows(enrich_c2_signif, df2) %>% 
  dplyr::select(pathway, group_in_pathway, FDR, genes, `k/K`) %>% 
  #Fix combined term with both C2 and C5 in it
  mutate(pathway = recode(pathway, 
                  "KEGG glutathione metabolism"="KEGG glutathione metabolism\nGOBP cellular oxidant detoxification",
                  "GOBP ribose phosphate biosynthetic process"="GOBP cellular carbohydrate metabolic process\nGOBP cellular ketone metabolic process\nGOBP ribose phosphate biosynthetic process")) %>% 
  filter(! pathway %in% c("GOBP cellular oxidant detoxification",
                          "GOBP cellular carbohydrate metabolic process",
                          "GOBP cellular ketone metabolic process"))



#### STRING network #####
map <- map_string(genes = genes.OI, score_threshold = 700)
p1 <- plot_string(map, enrichment = enrich_to_plot, layout = "mds",
                  colors = c("#CC79A7","#0072B2","#F0E442",
                             "#009E73","#56B4E9","#E69F00","grey80"),
                  text_size = 3, node_size = 2
                  ) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(fill = guide_legend(override.aes = list(size=13), title = "Parent term of enriched pathways"))
# p1

#### Save ####
ggsave(filename = "publication/archive/Fig4_enrichment.pdf", p1, width=12, height=6)
