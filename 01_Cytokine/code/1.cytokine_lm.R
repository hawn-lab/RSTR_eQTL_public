#Loading Required Packages and Setting Working Directory:
library(tidyverse)

#Date 
run_time <- gsub(":","_", Sys.Date())

#Reading and Filtering final eQTL Data:
eQTL_final <- read_csv("00_eQTL/results/2023-06-23_final_eQTL_definition.csv") %>% 
  filter(lead_SNP=="Y")

#cytokine gene
listOfCytokineGene <- c("TNF",
                        "IL6",
                        "IFNB1",
                        "IL1B"#,
                        #"CAMP",#not in rnaseq
                        #"DEFB1",#not in rnaseq
                        #"NOS2",#not in rnaseq
                        #"IRF3",#not used
                        #"NFKB1"#not used
                        )

#Gene expression data
attach("data_clean/eQTL_data_all.RData")

CytokineGene <- GENE %>%
  rownames_to_column("gene_hgnc_symbol") %>% 
  filter(gene_hgnc_symbol %in% listOfCytokineGene )

# Calculate delta
Cyto_delta <- CytokineGene %>% 
  pivot_longer(-gene_hgnc_symbol) %>% 
  separate(name, into=c("ptID","condition"), sep="_") %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(deltaTB = TB-MEDIA) %>% 
  select(-MEDIA, -TB) %>% 
  pivot_wider(names_from = ptID, values_from = deltaTB) %>% 
  column_to_rownames("gene_hgnc_symbol")

#Reformat explanatory variable to per patient
SNP_delta <- SNP %>% 
  select(contains("MEDIA")) %>% 
  rename_all(~gsub("_MEDIA","",.))
  
cov_n80 <- cov_n160  %>% 
  select(contains("MEDIA")) %>% 
  rename_all(~gsub("_MEDIA","",.))

# Initialize an empty data set
All_data_set <- c()

for(i in 1:nrow(Cyto_delta)){
  print(i)
  gene_name <- rownames(Cyto_delta)[i]
  # Loop through each row in eQTL
  for(j in 1:nrow(eQTL_final)) {
    # Get the SNP and gene names for the current row
    snps_name <- eQTL_final$snps[j]
    
    # Create a unique variable ID
    var_ID <- paste(snps_name, gene_name, sep = "_with_")
    
    # Extract the SNP and gene data
    snps_lm <- t(SNP_delta[snps_name,]) 
    colnames(snps_lm) <- "snp"
    gene_lm <- t(Cyto_delta[gene_name,])
    colnames(gene_lm)<-"expression"
    
    #Combine all model data
    dat <- full_join(rownames_to_column(as.data.frame(gene_lm)),
                     rownames_to_column(as.data.frame(snps_lm)),
                     by = join_by(rowname)) %>% 
      full_join(cov_n80 %>% t() %>% as.data.frame() %>% rownames_to_column(),
                by = join_by(rowname)) %>% 
      mutate(M0_KCVSEX = factor(M0_KCVSEX))
    
    # Define the index variable names
    index <- c("(Intercept)", snps_name, "SEX", "AGE")
    
    # Fit a linear model
    lmodel <- lm(expression ~ snp + M0_KCVSEX + M0_KCVAGE, data=dat)
    
    # Extract coefficient summary and create a data frame
    lmout <- summary(lmodel)$coefficients %>%
      data.frame() %>%
      dplyr::mutate(index = index) %>%
      dplyr::mutate(snps = snps_name) %>%
      dplyr::mutate(gene = gene_name) %>%
      dplyr::select(c("snps", "gene", "index", "Estimate", "t.value", "Pr...t..")) %>%
      dplyr::rename(c("p.value" = "Pr...t.."))
    
    # Append the selected row to the All_data_set data set
    All_data_set <- rbind(All_data_set, lmout["snp",])
  }
}

# Set the random seed for reproducibility
set.seed(13251246)

# Calculate the false discovery rate (FDR) adjusted p-values
All_data_set_fdr <- All_data_set %>% 
  mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
  arrange(FDR) %>% 
  #Add cis gene for SNP
  left_join(eQTL_final %>% rename(cis_gene=gene)) %>% 
  arrange(p.value)

# Save
dir.create("01_Cytokine/results/delta_model/", showWarnings = FALSE)

# Write the all data set to a CSV file
write_csv(All_data_set_fdr, 
          file=paste0("01_Cytokine/results/delta_model/", run_time,
                      "_cytokine_interaction_co.csv"))

