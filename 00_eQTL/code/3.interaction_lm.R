#Loading Required Packages and Setting Working Directory:
library(tidyverse)

#Load data
# Created in 00_eQTL/data_cleaning.R
load("data_clean/eQTL_data_all.RData")

#Date 
run_time <- gsub(":","_", Sys.Date())

# Check data order
identical(colnames(GENE), colnames(SNP))
identical(colnames(GENE), colnames(cov_n160))

#Load eQTL results
eQTL <- read_csv("00_eQTL/results/Matrix_eQTL/2023-06-14_eqtl_TB_co.csv",
                 col_types = cols(
                   snps = col_character(),
                   gene = col_character(),
                   statistic = col_double(),
                   pvalue = col_double(),
                   FDR = col_double(),
                   beta = col_double()))%>%
  dplyr::filter(FDR < 0.01) %>%
  select(snps, gene)

# Initialize an empty data set
All_data_set <- c()

# Loop through each row in eQTL
for(i in 1:nrow(eQTL)) {
  print(i)
  # Get the SNP and gene names for the current row
  snps_name <- eQTL$snps[i]
  gene_name <- eQTL$gene[i]
  
  # Create a unique variable ID
  var_ID <- paste(snps_name, gene_name, sep = "_with_")
  
  # Extract the SNP and gene data
  snps_lm <- t(SNP[snps_name,]) 
  colnames(snps_lm) <- "snp"
  gene_lm <- t(GENE[gene_name,])
  colnames(gene_lm)<-"expression"
  
  #Combine all model data
  dat <- full_join(rownames_to_column(as.data.frame(gene_lm)),
                   rownames_to_column(as.data.frame(snps_lm)),
                   by = join_by(rowname)) %>% 
    separate(rowname, into = c("RSID","stim"), sep="_", remove = FALSE) %>% 
    full_join(cov_n160 %>% t() %>% as.data.frame() %>% rownames_to_column,
              by = join_by(rowname)) %>% 
    mutate(M0_KCVSEX = factor(M0_KCVSEX))
  
  # Define the index variable names
  index <- c("(Intercept)", snps_name, "stim_status", "SEX", "AGE", paste(snps_name, "stim_status", sep = ":"))
  
  # Fit a linear model
  lmodel <- lm(expression ~ snp*stim + M0_KCVSEX + M0_KCVAGE, data=dat)
  
  # Extract coefficient summary and create a data frame
  lmout <- summary(lmodel)$coefficients %>%
    data.frame() %>%
    dplyr::mutate(index = index) %>%
    dplyr::mutate(snps = snps_name) %>%
    dplyr::mutate(gene = gene_name) %>%
    dplyr::select(c("snps", "gene", "index", "Estimate", "t.value", "Pr...t..")) %>%
    dplyr::rename(c("p.value" = "Pr...t.."))
  
  # Append the selected row to the All_data_set data set
  All_data_set <- rbind(All_data_set, lmout[6,])
}

# Set the random seed for reproducibility
set.seed(13251246)

# Calculate the false discovery rate (FDR) adjusted p-values
All_data_set_fdr <- All_data_set %>% 
  mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
  arrange(FDR)

# Save
dir.create("00_eQTL/results/interaction_model/", showWarnings = FALSE)

# Write the all data set to a CSV file
write_csv(All_data_set_fdr, 
          file=paste0("00_eQTL/results/interaction_model/", run_time,
                      "_eQTL_interaction_co.csv"))

