#Loading Required Packages and Setting Working Directory:
library(tidyverse)

#### SNP DATA ####
snp_original <- read_tsv("data/snp/snp_unfiltered_MAF_HWE_10e_6_n80_filtered_by_annot_201118.txt.gz")

#Rename to FULLIDNO
RSID <- colnames(snp_original)[-1]
ID_key <- read_csv("data/metadata/2023.03.13RSTR_Hawn_metadata.csv") %>% 
  select(FULLIDNO,RS_SUB_ACCESSION_NO) %>% 
  drop_na() %>% 
  distinct()
FULLID <- ID_key %>% 
  filter(RS_SUB_ACCESSION_NO %in% RSID) %>% 
  arrange(factor(RS_SUB_ACCESSION_NO, levels = RSID))
#Check order
identical(FULLID$RS_SUB_ACCESSION_NO, RSID)

#Rename
snp_rename <- snp_original
  colnames(snp_rename) <- c("rsID", FULLID$FULLIDNO)
#Reorder
FULLID_ord <- colnames(snp_rename)[-1] %>% sort()
snp_rename_ord <- snp_rename[,c("rsID",FULLID_ord)]

#copy for use with media and tb rnaseq samples
snp_tb <- snp_rename_ord
  snp_tb_names <- paste(colnames(snp_tb)[-1],"TB",sep="_")
  colnames(snp_tb) <- c("rsID",snp_tb_names)

snp_media <- snp_rename_ord
  snp_media_names <- paste(colnames(snp_media)[-1],"MEDIA",sep="_")
  colnames(snp_media) <- c("rsID",snp_media_names)

all_ID  <- str_sort(c(snp_media_names, snp_tb_names))

# Combine SNP data
SNP <- snp_media %>%
  left_join(snp_tb, by = "rsID") %>%
  dplyr::select(rsID, all_of(all_ID))

#### RNASEQ DATA ####
# Transform gene expression data
gene_original <- read_tsv("data/rnaseq/GE.160_after_QC_voom_n80_201121.txt",
                 col_types = cols(
                   .default = col_double(),
                   gene_hgnc_symbol = col_character()))

# Rename to FULLID
RSID2 <- vapply(strsplit(colnames(gene_original)[-1], "_", fixed = TRUE), "[", "", 1) %>% 
  unique()
FULLID2 <- ID_key %>% 
  filter(RS_SUB_ACCESSION_NO %in% RSID2)

#Rename
GENE <- gene_original %>% 
  pivot_longer(-gene_hgnc_symbol) %>% 
  separate(name, into=c("RS_SUB_ACCESSION_NO","mtb"), sep="_") %>% 
  left_join(FULLID2) %>% 
  mutate(name=paste(FULLIDNO,mtb,sep="_")) %>% 
  select(-RS_SUB_ACCESSION_NO, -FULLIDNO, -mtb) %>% 
  arrange(factor(name, levels=all_ID)) %>% 
  pivot_wider()

identical(colnames(GENE)[-1],all_ID)

#### Variables for model ####
# Read covariates from the file
covariates <- read_csv("data/metadata/2023.03.13RSTR_Hawn_metadata.csv") %>% 
  distinct(FULLIDNO, M0_KCVSEX, M0_KCVAGE) %>%  #KCHCA_AGE_YR_CURRENT
  #recode sex
  mutate(M0_KCVSEX = case_when(M0_KCVSEX=="M"~0,
                               M0_KCVSEX=="F"~1)) %>% 
  filter(FULLIDNO %in% FULLID$FULLIDNO) %>% 
  arrange(factor(FULLIDNO, levels=FULLID_ord)) %>% 
  column_to_rownames("FULLIDNO") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("index")

#Make media and tb covariates
cov_tb <- covariates
  cov_tb_names <- paste(colnames(cov_tb)[-1],"TB",sep="_")
  colnames(cov_tb) <- c("index",cov_tb_names)

cov_media <- covariates
  cov_media_names <- paste(colnames(cov_media)[-1],"MEDIA",sep="_")
  colnames(cov_media) <- c("index",cov_media_names)

# Combine SNP data
cov_n160 <- cov_media %>%
  left_join(cov_tb, by = "index") %>%
  dplyr::select(index, all_of(all_ID))

identical(colnames(cov_n160)[-1], all_ID)

#### Save ####
dir.create("data_clean/", showWarnings = FALSE)
SNP <- SNP %>% column_to_rownames("rsID")
GENE <- GENE %>% column_to_rownames("gene_hgnc_symbol")
cov_n160 <- cov_n160 %>% column_to_rownames("index")
save(SNP, GENE, cov_n160, file = "data_clean/eQTL_data_all.RData")

#Split media and tb
GENE_MEDIA <- GENE %>% select(contains("MEDIA"))
GENE_TB <- GENE %>% select(contains("TB"))

SNP_MEDIA <- SNP %>% select(contains("MEDIA"))
SNP_TB <- SNP %>% select(contains("TB"))

COV_MEDIA <- cov_n160 %>% select(contains("MEDIA"))
COV_TB <- cov_n160 %>% select(contains("TB"))

save(SNP_MEDIA, GENE_MEDIA, COV_MEDIA,
     SNP_TB, GENE_TB, COV_TB,
     file = "data_clean/eQTL_data_by_stim.RData")

# write.table(SNP, sep="\t", 
#             row.names=FALSE, col.names=TRUE,
#             file=gzfile("data_clean/SNP_data.txt.gz"))
# write_tsv(GENE, file="data_clean/GENE_data.txt")
# write_tsv(GENE_MEDIA, file="data_clean/GENE_media_data.txt")
# write_tsv(GENE_TB, file="data_clean/GENE_tb_data.txt")
