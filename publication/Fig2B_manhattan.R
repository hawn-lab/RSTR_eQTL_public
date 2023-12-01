library(tidyverse)
library(ggrepel)

#### Data ###
# SNP position
POS <- read_tsv("data/snp/snpsloc_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt") %>% 
  rename(snps=id)

# Mtb eQTL
eQTL_mtb <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv") %>% 
  #Add genome position
  left_join(POS) %>% 
  select(snps, gene, chr, pos, FDR) %>% 
  mutate(mtb = ifelse(FDR < 0.01, "Y",NA))

#Media eQTL
eQTL_media <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv") %>% 
  filter(FDR < 0.01) %>% 
  select(snps, gene) %>% 
  mutate(media="Y")

# Interaction eQTL
eQTL_inter <- read_csv("00_eQTL/results/interaction_model/2023-06-21_eQTL_interaction_co.csv") %>% 
  filter(FDR < 0.1) %>% 
  select(snps, gene) %>% 
  mutate(inter="Y")

#### Format data ####
#Identify signif eQTL
dat <- eQTL_mtb %>% 
  left_join(eQTL_inter) %>% 
  left_join(eQTL_media) %>% 
  #Create alternating chr color groups
  mutate(color = case_when(mtb == "Y" & inter=="Y" & is.na(media) ~ "Significant cis eQTL",
                           chr %in% seq(1,22,2) ~ "Odd chr",
                           chr %in% seq(2,22,2) ~ "Even chr"))

#Modified from ggman::ggman( )
#Add initial index in order
dat_index <- dat[order(dat$pos), ]
dat_index <- dat_index[gtools::mixedorder(dat_index$chr), ]
dat_index$index <- 1:nrow(dat_index)

relpos <- function(x, minbp, maxbp, nrows, startingpoint) {
  actual.distance <- (x - minbp)
  relative.distance <- (actual.distance * nrows)/maxbp
  return(relative.distance + startingpoint) 
}

#Count SNP in chr
chrtable <- data.frame(table(dat_index$chr))
dat_index$chr <- factor(dat_index$chr, levels = chrtable$Var1)
dfm.split <- split(dat_index, dat_index$chr)
startingpoint = 1
dfm.list <- lapply(dfm.split, function(x) {
  minbp <- as.numeric(min(x$pos))
  maxbp <- as.numeric(max(x$pos))
  nrows <- as.numeric(nrow(x))
  x$index <- relpos(x$pos, minbp, maxbp, nrows, startingpoint)
  startingpoint <<- max(x$index) + 1
  return(x)
})

dat_index_rel <- do.call(rbind, dfm.list) %>% 
  arrange(color)

#Create label group
gene_labels <- dat_index_rel %>% 
  mutate(label = case_when(mtb == "Y" & inter=="Y" & is.na(media) ~ gene)) %>% 
  drop_na(label) %>% 
  distinct(snps, label, FDR, index) %>% 
  group_by(label) %>% 
  summarise(index = min(index),
            FDR = min(FDR)) %>% 
  ungroup() %>% 
  distinct(label, FDR, index)

#Create Chr labels
xbreaks <- dat_index_rel %>% 
  group_by(chr) %>% 
  summarise(br = mean(index)) %>% 
  ungroup()

#### Plot ####

p1 <-dat_index_rel %>% 
  
  ggplot(aes(x=index, y=-log10(FDR)))  +
  #label genes
  geom_text_repel(data=gene_labels, aes(label = label), color="red",
                  force=1, point.padding=unit(1,'lines'),
                  vjust=1, direction='y', 
                  nudge_x=0,nudge_y = 20-(-log10(gene_labels$FDR)),
                  segment.size=0.2, segment.color="grey80") +
  #Add points
  geom_point(aes(color=color, size=color))+
  scale_x_continuous(breaks = xbreaks$br, labels = xbreaks$chr, expand = c(0, 0)) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.01), lty="dashed") +
  scale_color_manual(values=c("grey60","grey20","red")) +
  scale_size_manual(values=c(1,1,3)) +
  labs(x="Chromosome")+
  theme(legend.position = "none")

ggsave("publication/Fig2B_manhattan.png", p1, width=11, height=4)
ggsave("publication/Fig2B_manhattan.pdf", p1, width=11, height=4)
