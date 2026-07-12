# SCZ mQTL eQTL Colocalisation Pipeline
# Script: eQTL Preprocessing and Colocalisation
# Author: Ruthushree P Basappa
# Institution: University of Galway
# Supervisor: Prof. Derek Morris
# Date: June 2026

library(data.table)
library(dplyr)
library(here)

here::i_am("01b_data_preprocessing_eqtl/eqtl_R_analysis.R")

# Step 1: Extract SNPs within GWAS loci from full O'Brien eQTL dataset
# (all_eqtls_gene.txt obtained from figshare: https://doi.org/10.6084/m9.figshare.6881825)

loci_snps <- fread(here("data", "loci_snp_list.csv"))
my_snps <- loci_snps$SNP

eqtl_loci <- fread(here("data", "all_eqtls_gene.txt"))[variant_id %in% my_snps]
cat("eQTL rows extracted:", nrow(eqtl_loci), "\n")

fwrite(eqtl_loci, here("data", "eqtl_loci_only.csv"))

# Step 2: Rename columns to match COLOC-reporter conventions

eqtl_loci <- eqtl_loci %>% rename(
  GeneID       = gene_id,
  SNP          = variant_id,
  TSS_distance = tss_distance,
  MA_samples   = ma_samples,
  MA_count     = ma_count,
  MAF          = maf,
  P_eqtl       = pval_nominal,
  BETA_eqtl    = slope,
  SE_eqtl      = slope_se
)

# Step 3: Add CHR and POS from reference SNP list

eqtl_loci <- inner_join(
  eqtl_loci,
  loci_snps %>% dplyr::select(CHR, POS, SNP, A1, A2) %>%
    rename(A1_ref = A1, A2_ref = A2),
  by = "SNP"
)
cat("After adding CHR and POS:", nrow(eqtl_loci), "\n")

# Step 4: QC and sample size / MAF assignment

eqtl_loci <- eqtl_loci %>%
  filter(is.finite(SE_eqtl), SE_eqtl > 0) %>%
  mutate(N_eqtl = 120)

gwas <- fread(here("data", "gwas_qc.csv"))
gwas_maf <- gwas %>% dplyr::select(SNP, MAF)
rm(gwas); gc()

eqtl_loci <- eqtl_loci %>%
  rename(MAF_eqtl = MAF) %>%
  left_join(gwas_maf, by = "SNP") %>%
  mutate(MAF = ifelse(is.na(MAF), MAF_eqtl, MAF))

fwrite(eqtl_loci, here("data", "eqtl_preprocessed.csv"))
cat("Preprocessed eQTL saved\n")

# Step 5: Locus-level colocalisation was run using COLOC-reporter across all
# 201 independent GWAS loci (see 02b_colocalisation_eqtl/ for the bash script
# and GWAS_samples_eqtl.txt / set.regions.eqtl.txt configuration files).
#
# Loci passing PP4 >= 0.8 at the locus level were then carried forward for
# per-gene colocalisation following Dobbyn et al. (2018), run independently
# for each gene with available eQTL data within the locus window.
