# Prepare GWAS summary statistics for SMR analysis
# Converts rsID format to chr:bp format to match FB_Brain_2 BESD file
# Author: Ruthushree P Basappa
# Date: May 2026
library(data.table)
library(here)
here::i_am("04_SMR_analysis/prepare_gwas_smr.R")
dir.create(here("results", "smr"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("references", "g1000_eur_chrbp"), showWarnings = FALSE, recursive = TRUE)
gwas_qc <- fread(here("data", "gwas_qc.csv"))
gwas_qc[, SNP_chrbp := paste0(CHR, ":", POS)]
gwas_smr_chrbp <- gwas_qc[, .(
  SNP = SNP_chrbp,
  A1 = EA,
  A2 = OA,
  freq = MAF,
  b = BETA_gwas,
  se = SE_gwas,
  p = P_gwas,
  N = NEFF
)]
fwrite(gwas_smr_chrbp, here("data", "gwas_smr_chrbp.txt"), sep = "\t")
cat("GWAS SMR file created:", nrow(gwas_smr_chrbp), "SNPs\n")
source_prefix <- here("references", "g1000_eur", "g1000_eur")
target_prefix <- here("references", "g1000_eur_chrbp", "g1000_eur")
system(sprintf("cp '%s.fam' '%s.fam'", source_prefix, target_prefix))
system(sprintf("cp '%s.bed' '%s.bed'", source_prefix, target_prefix))
system(sprintf(
  "awk '{print $1\"\\t\"$1\":\"$4\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}' '%s.bim' > '%s.bim'",
  source_prefix, target_prefix
))
cat("LD reference converted to chr:bp format\n")
