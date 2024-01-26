# Set required libraries
library(tidyverse)    # Collection of data manipulation and visualization packages

#### eQTL data ####
#Final TB induced eQTL
# filtered_eQTL <- read_csv("00_data/filtered_eQTL_Oct21_2021.csv")
eQTL <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv")

#TB eQTLs
TB_eQTL <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv",
                    col_types = cols(
                      snps = col_character(),
                      gene = col_character(),
                      statistic = col_double(),
                      pvalue = col_double(),
                      FDR = col_double(),
                      beta = col_double()
                    )) %>%
  dplyr::filter(FDR < 0.01) %>%
  arrange(pvalue) %>%
  distinct(snps, .keep_all = TRUE)

#Media eQTLs
MEDIA_eQTL <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_MEDIA_co.csv",
                       col_types = cols(
                         snps = col_character(),
                         gene = col_character(),
                         statistic = col_double(),
                         pvalue = col_double(),
                         FDR = col_double(),
                         beta = col_double()
                       )) %>%
  dplyr::filter(FDR < 0.01) %>%
  arrange(pvalue) %>%
  distinct(snps, .keep_all = TRUE)

#### Genotype data ###
#Merging SNP and Gene Data:
snp_original <- read_tsv("data/snp/snp_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt.gz")

#Create separate media data
snp_media <- snp_original

snp_media_names <- paste(substr(colnames(snp_media[,2:ncol(snp_media)]),
                                regexpr('R', colnames(snp_media[,2:ncol(snp_media)])), 
                                regexpr('R', colnames(snp_media[,2:ncol(snp_media)]))+7),
                         "MEDIA",sep = "_")

colnames(snp_media) <- c("rsID",snp_media_names)

#Create separate TB data
snp_tb <- snp_original
snp_tb_names <- paste(substr(colnames(snp_tb[,2:ncol(snp_tb)]),
                             regexpr('R', colnames(snp_tb[,2:ncol(snp_tb)])), 
                             regexpr('R', colnames(snp_tb[,2:ncol(snp_tb)]))+7),
                      "TB",sep = "_")
colnames(snp_tb) <- c("rsID",snp_tb_names)

# Combine media and TB data
##Order by RSID
all_RSID  <- str_sort(c(snp_media_names, snp_tb_names))

SNP <- snp_media %>%
  left_join(snp_tb, by = "rsID") %>%
  data.frame() %>%
  column_to_rownames(var = "rsID") %>%
  dplyr::select(all_of(all_RSID))


#### RNAseq data ####
# Read gene expression data
# Transform gene expression data
GENE <- read_tsv("data/RNAseq/GE.160_after_QC_voom_n80_201121.txt",
                 col_types = cols(
                   .default = col_double(),
                   gene_hgnc_symbol = col_character()
                 )) %>%
  data.frame() %>%
  column_to_rownames(var = "gene_hgnc_symbol") %>%
  dplyr::select(all_of(all_RSID))

#### Covariate data ####
# Create 01 stimulus status matrix
stim <-  data.frame(RSID = all_RSID) %>%
  mutate(stim_status = ifelse(grepl("_MEDIA", RSID), 0, 1)) %>%
  mutate(stim_status = as.numeric(stim_status)) %>% 
  column_to_rownames(var = "RSID") %>% 
  as.matrix()

# Read covariate data
# Read covariates from the file
cov <- read_csv("data/metadata/2023.03.13RSTR_Hawn_metadata.csv") %>% 
  distinct(RS_SUB_ACCESSION_NO, M0_KCVSEX, M0_KCVAGE) %>% 
  #recode sex
  mutate(M0_KCVSEX = case_when(M0_KCVSEX=="M"~0,
                               M0_KCVSEX=="F"~1)) %>% 
  #Keep samples in SNP data
  filter(RS_SUB_ACCESSION_NO %in% colnames(snp_original))

cov_n160 <- cov %>% 
  mutate(RS_SUB_ACCESSION_NO=paste(RS_SUB_ACCESSION_NO,"MEDIA",sep="_")) %>% 
  bind_rows(cov %>% 
              mutate(RS_SUB_ACCESSION_NO=paste(RS_SUB_ACCESSION_NO,"TB",sep="_"))) %>% 
  arrange(RS_SUB_ACCESSION_NO) %>% 
  #matrix format
  column_to_rownames("RS_SUB_ACCESSION_NO") %>% 
  as.matrix()

#### Run models ####
# Initialize an empty data set
All_data_set <- c()

# Loop through each row in eQTL
for(i in 1:nrow(eQTL)) {
  # Get the SNP and gene names for the current row
  snps_name <- eQTL$snps[i]
  gene_name <- eQTL$gene[i]
  
  # Create a unique variable ID
  var_ID <- paste(snps_name, gene_name, sep = "_with_")
  
  # Extract the SNP and gene data
  snps_lm <- t(SNP[snps_name,]) 
  gene_lm <- t(GENE[gene_name,])
  
  #Check orders
  snp_ord <- identical(rownames(snps_lm), rownames(cov_n160))
  gene_ord <- identical(rownames(gene_lm), rownames(cov_n160))
  if(any(!snp_ord, !gene_ord)){stop("Data order mismatch")}
  
  # Define the index variable names
  index <- c("(Intercept)", snps_name, "stim_status", "SEX", "AGE", paste(snps_name, "stim_status", sep = ":"))
  
  # Fit a linear model
  lmodel_dominant <- lm(gene_lm ~ I(snps_lm < 2) * stim + cov_n160)
  lmodel_recessive <- lm(gene_lm ~ I(snps_lm == 2) * stim + cov_n160)
  lmodel_additive <- lm(gene_lm ~ snps_lm  * stim + cov_n160)
  
  #Calculate AIC
  AIC_dominant <- AIC(lmodel_dominant)
  AIC_recessive <- AIC(lmodel_recessive)
  AIC_additive <- AIC(lmodel_additive)
  
  #Save
  All_data_set <- data.frame(gene = c(gene_name),
                             snp = c(snps_name),
                             additive = c(AIC_additive),
                             recessive = c(AIC_recessive),
                             dominant = c(AIC_dominant)) %>% 
    bind_rows(All_data_set)
}

#Calculate AIC differences
All_data_set <- All_data_set %>% 
  mutate(AIC_dominant = additive-dominant,
         AIC_recessive = additive-recessive)

# Write the final data set to a CSV file
write_csv(All_data_set, "00_eQTL/results/model_AIC_dominant_recessive.csv")

