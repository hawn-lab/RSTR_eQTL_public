library(tidyverse)
library(genetics)
library(snpStats)
library(ggrepel)
library(patchwork)
# gene.OI <- c("SIRPB1","PRKAG2","NDUFAF4")

#### Data ####
attach("data_clean/eQTL_data_all.RData")
eqtl <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv")
eqtl_lm <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv")
eqtl_key <- read_csv("data/snp/eQTL_rsID_key.csv")
snp_pos <- read_tsv("data/snp/snpsloc_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt")
snp_key <- read_tsv("data/snp/snp.annot_megaEX_unfiltered_n_80_filtered_by_annot_201118.txt")

#### Calculate LD and plot ####

plot.ls <- list()

for(gene.OI in c("SIRPB1","PRKAG2","NDUFAF4")){
  print("Calculating LD")
  eqtl_filter <- eqtl %>% filter(gene == gene.OI)
  
  #Get all eQTL tested for gene of interest
  eqtl_lm_filter <- eqtl_lm %>% filter(gene == gene.OI) %>% 
    left_join(eqtl_key)
  
  #Get SNP positions
  snp_pos_filter <- snp_pos %>% filter(id %in% eqtl_lm_filter$snps)
  
  #### Calculate LD ####
  # Filter to SNP of interest and convert to alleles
  snp_key_filter <- snp_key %>% filter(rsID %in% snp_pos_filter$id)
  
  snp_filter <- SNP %>% 
    rownames_to_column("snpID") %>% 
    filter(snpID %in% snp_pos_filter$id) %>% 
    pivot_longer(-snpID) %>%
    left_join(snp_key_filter, by=c("snpID"="rsID")) %>%
    mutate(value = case_when(value==0 ~ paste(ref,ref,sep="/"),
                             value==1 ~ paste(ref,alt,sep="/"),
                             value==2 ~ paste(alt,alt,sep="/"))) %>%
    dplyr::select(-ref,-alt) %>%
    pivot_wider()
  
  #Format values and order SNPs in chr
  snp_ord <- snp_pos_filter %>% arrange(pos) %>% pull(id)
  snp_filter_ord <- snp_filter %>%
    mutate(snpID = factor(snpID, levels=snp_ord)) %>%
    arrange(snpID)

  # List all SNP pairs
  snps.to.LD <- t(combn(snp_filter_ord$snpID, m=2)) %>%
    as.data.frame() %>%
    #Remove duplicates
    distinct(V1,V2) %>%
    filter(V1!=V2)

  #Calculate LD
  LD.df <- snps.to.LD %>%
    mutate(R2="")

  for(i in 1:nrow(LD.df)){
    print(i)
    #Get SNP pair for i index
    SNP.pair <- snps.to.LD[i, ] %>% unlist(use.names = FALSE)

    #Filter allele data to SNP pairs
    snp_filter.sub <- snp_filter %>%
      filter(snpID %in% SNP.pair) %>%
      column_to_rownames("snpID") %>%
      t()

    #Calculate LD R^2
    #Only run if > 1 alleles in bpth SNPs
    no.snp1 <- unique(snp_filter.sub[,1][!is.na(snp_filter.sub[,1])]) %>%
      length()
    no.snp2 <- unique(snp_filter.sub[,2][!is.na(snp_filter.sub[,2])]) %>%
      length()

    if(no.snp1 > 1 & no.snp2 > 1){
      LD.df[i, "R2"] <- genetics::LD(as.genotype(snp_filter.sub[,1]),
                                     as.genotype(snp_filter.sub[,2]))[["R^2"]]
    } else{
      LD.df[i, "R2"] <- NA
    }
  }

  LD.df.format <- LD.df %>%
    mutate(R2=as.numeric(R2)) %>%
    mutate(V1 = factor(V1, levels=snp_ord),
           V2 = factor(V2, levels=snp_ord))

  save(LD.df.format, file=paste0("publication/data/",gene.OI,".LD.RData"))

  # load(paste0("publication/data/",gene.OI,".LD.RData"))
  
  print("Plotting")
  #### Heatmaps ####
  LD3 <- LD.df.format %>% 
    
    ggplot(aes(x = V1, y = V2, fill = R2)) +
    geom_raster() +
    scale_fill_gradient2(low="white",mid = "white",high="red") +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          aspect.ratio = 1) +
    coord_flip()
  # LD3
  
  #### eQTL significance ####
  eqtl_pos <- eqtl_lm_filter %>%
    mutate(snps = factor(snps, levels=snp_ord)) %>% 
    left_join(eqtl_filter %>% dplyr::select(snps,gene,lead_SNP,rsID)) %>% 
    left_join(snp_pos_filter, by=c("snps"="id"))
  
  LD1 <- eqtl_pos %>% 
    mutate(color = case_when(snps %in% eqtl$snps ~ "eQTL",
                             TRUE ~ "Not significant"),
           label = case_when(lead_SNP =="Y" ~ rsID)) %>% 
    
    ggplot(aes(x=pos, y=-log10(FDR))) +
    geom_point(alpha=0.5, shape=16) +
    geom_text_repel(aes(label = label), direction = "x", show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.01), lty="dashed") +
    theme_classic() +
    labs(x = "", color="") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  # LD1
  
  #### Connecting lines ####
  span <- max(eqtl_pos$pos)-min(eqtl_pos$pos)
  interval <- span/(length(eqtl_pos$pos)-1)
  
  LD2 <- data.frame(top = sort(eqtl_pos$pos),
                    bot = seq(min(eqtl_pos$pos), max(eqtl_pos$pos), interval)) %>% 
    
    ggplot() +
    geom_segment(aes(x=bot, xend=top), y=0, yend=1) +
    theme_void()

  plot.ls[[gene.OI]] <- LD1+LD2+LD3 +
    plot_layout(ncol=1, heights = c(0.3,0.1,sqrt(0.5)))

}

#### Save ####

ggsave(filename = "publication/archive/Fig3_LD.heatmap.pdf", plot.ls[["SIRPB1"]], 
       width=6, height=6)

ggsave(filename = "publication/archive/FigS2_LD.heatmap1.pdf", 
       plot.ls[["PRKAG2"]]+plot_annotation(title="a.  PRKAG2"), 
       width=6, height=6)
ggsave(filename = "publication/archive/FigS2_LD.heatmap2.pdf", 
       plot.ls[["NDUFAF4"]]+plot_annotation(title="b.  NDUFAF4"),  
       width=6, height=6)
