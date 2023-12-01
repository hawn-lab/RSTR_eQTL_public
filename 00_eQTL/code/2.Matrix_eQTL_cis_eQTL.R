# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r")
library(MatrixEQTL)
library(tidyverse)

set.seed(497303827)

#Date 
run_time <- gsub(":","_", Sys.Date())

#Make folder
dir.create("00_eQTL/results/Matrix_eQTL/",
           recursive = TRUE, showWarnings = FALSE)

#Attach clean data
attach("data_clean/eQTL_data_by_stim.RData")

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
#useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# SNP location file name
snps_location_file_name = "data/snp/snpsloc_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt"

# Gene location file name
gene_location_file_name = "data/rnaseq/gene_location_final_201118.txt"

# Annotation File
SNP_annot_file_name = "data/snp/snp.annot_megaEX_unfiltered_n_80_filtered_by_annot_201118.txt"

# Output file name
output_file_name_cis = NULL
output_file_name_tra = NULL

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()
# errorCovariance = read.table("Sample_Data/errorCovariance.txt")

# Distance for local gene-SNP pairs
cisDist = 1e6

#Loop for MEDIA and TB analysis
for(stim in c("MEDIA","TB")){
  print(stim)
  
  if(stim=="MEDIA"){ file_name <- "_eqtl_MEDIA_co" }
  if(stim=="TB"){ file_name <- "_eqtl_TB_co" }
  
  ## Load genotype data
  snps <- SlicedData$new()
  snps$CreateFromMatrix(as.matrix(get(paste0("SNP_",stim))))

  ## Load covariates
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(get(paste0("COV_",stim))))
  
  ## Load gene expression data
  gene <- SlicedData$new()
  gene$CreateFromMatrix(as.matrix(get(paste0("GENE_",stim))))
  
  ## Check data order
  check1 <- identical(gsub("\"", "", snps$columnNames), gene$columnNames)
  check2 <- identical(gsub("\"", "", snps$columnNames), cvrt$columnNames)
  if(any(!check1, !check2)){stop("Data are not in the same column order.")}
  
  ## Run the analysis
  snpspos = read_tsv(snps_location_file_name,
                     col_types = cols(
                       id = col_character(),
                       chr = col_double(),
                       pos = col_double()
                     )) %>%
    as.data.frame()
  
  genepos <- read_tsv(gene_location_file_name,
                     col_types = cols(
                       gene_hgnc_symbol = col_character(),
                       chromosome_name = col_character(),
                       start_position = col_double(),
                       end_position = col_double()
                     )) %>%
    as.data.frame() %>% 
    filter(gene_hgnc_symbol %in% rownames(get(paste0("GENE_",stim))))
  
  me <- Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = "",
    pvOutputThreshold = 0,
    useModel = modelLINEAR, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = "",
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE)
  
  unlink(output_file_name_cis)
  
  ## Results:
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')

  ## Plot the Q-Q plot of local and distant p-values
  dir.create("00_eQTL/figs", showWarnings = FALSE)
  png(paste0("00_eQTL/figs/QQplot_Matrix",file_name,".png"), width=400,height = 400)
  plot(me)
  dev.off()
  
  #Save Results
  eqtl.cis <- me$cis$eqtls %>%
    data.frame()
  
  #save file
  filename.csv = paste(run_time,file_name,".csv",sep="")
  write_csv(eqtl.cis, paste0("00_eQTL/results/Matrix_eQTL/", filename.csv))
}
