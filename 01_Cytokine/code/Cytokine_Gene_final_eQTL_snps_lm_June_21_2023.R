#call library
library(tidyverse)

#set work directory
setwd("/Users/hyejeonghong/Library/Mobile Documents/com~apple~CloudDocs/03_Hyejeong_DATA/05_Code/03_Cytokine/")

####make snps data set####
#call snp data
snp <- read_tsv("data/snp_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt") %>%
  column_to_rownames("rsID") %>%
  filter()
#call rsID data
rsID_from_all_variants_file <- readxl::read_xlsx("data/Toms_request_all_variants_1MB_flanking_073120.xlsx") %>%
  dplyr::rename(c("snp" = "variant.id")) %>%
  dplyr::select(c("snp", "rsID"))

#call eqtl list
eqtls <- read_csv("data/results_eqtl_with_annot_and_rsID_Jan_16_2022.csv") %>%
  dplyr::select(SNPS) 

snp_RSID <- colnames(snp)



####make covariate data set####
cvrt <- read_tsv("data/covariates_SEX_AGE_201121.txt") %>%
  data.frame() %>%
  column_to_rownames(var = "index") %>%
  select(snp_RSID) %>%
  t()

#cytokine gene
listOfCytokineGene <- data.frame(gene_hgnc_symbol = c("TNF",
                                                      "IL6",
                                                      "IFNB1",
                                                      "IL1B",
                                                      "CAMP",
                                                      "DEFB1",
                                                      "NOS2",
                                                      "IRF3",
                                                      "NFKB1")) 

Gene_delta <- read_tsv("data/GE.delta_after_QC_voom_0_1_n80_201121.txt",
                       col_types = cols(
                         .default = col_double(),
                         gene_hgnc_symbol = col_character()
                       )) %>%
  column_to_rownames(var = "gene_hgnc_symbol") %>%
  select(snp_RSID)




# Initialize an empty data set
All_data_set <- c()

for(n in 1:nrow(listOfCytokineGene)){
  gene_name <- listOfCytokineGene$gene_hgnc_symbol[n]
  # Loop through each row in eqtls
  for(i in 1:nrow(eqtls)) {
    # Get the SNP and gene names for the current row
    snps_name <- eqtls$SNPS[i]
    
    
    # Create a unique variable ID
    var_ID <- paste(snps_name, gene_name, sep = "_with_")
    
    # Extract the SNP and gene data
    snps_lm <- t(snp [snps_name,]) 
    gene_lm <- t(Gene_delta[gene_name,])
    
    # Define the index variable names
    index <- c("(Intercept)", snps_name, "SEX", "AGE")
    
    # Fit a linear model
    lmodel <- lm(gene_lm ~ snps_lm + cvrt, na.action=) #lm(gene_lm ~ snps_lm + cvrt)
    
    # Extract coefficient summary and create a data frame
    lmout <- summary(lmodel)$coefficients %>%
      data.frame() %>%
      #  dplyr::mutate(index = index) %>%
      dplyr::mutate(snps = snps_name) %>%
      dplyr::mutate(gene = gene_name) %>%
      #     dplyr::select(c("snps", "gene", "index", "Estimate", "t.value", "Pr...t..")) %>%
      dplyr::rename(c("p.value" = "Pr...t..")) %>%
      dplyr::tibble()
    
    # Append the selected row to the All_data_set data set
    All_data_set <- rbind(All_data_set, lmout[2,])
  }
  
  
  
}

# Set the random seed for reproducibility
set.seed(13251246)

# Convert p-values to numeric
p_value <- All_data_set$p.value %>%
  as.numeric()

# Calculate the false discovery rate (FDR) adjusted p-values
All_data_set$FDR <- p.adjust(p_value, method = "BH")

All_data_set <- All_data_set %>%
  arrange(gene)%>%
  filter(p.value < 0.05)

# Create the final data set by arranging rows based on FDR, filtering for FDR < 0.1, and creating a tibble
data_set_final <- All_data_set %>%
  arrange(FDR) %>%
  dplyr::filter(FDR < 0.1) %>%
  tibble()

# Create a separate data set containing all rows, arranged by FDR, as a tibble
data_set_final_all <- All_data_set %>%
  arrange(FDR) %>%
  tibble()

# Keep only distinct rows based on the 'gene' column in the final data set
data_set_final %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Write the final data set to a CSV file
write_csv(data_set_final, "results/Cytokine_Gene_final_eQTL_snps_lm_FDR_0_1_cov_June_21_2023.csv")

# Write the final data set to a tab-delimited text file
write.table(data_set_final, "results/Cytokine_Gene_final_eQTL_snps_lm_FDR_0_1_cov_June_21_2023.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# Write the all data set to a CSV file
write_csv(data_set_final_all, "results/Cytokine_Gene_final_eQTL_snps_lm_June_21_2023.csv")

# Write the all data set to a tab-delimited text file
write.table(data_set_final_all, "results/Cytokine_Gene_final_eQTL_snps_lm_cov_June_21_2023.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

