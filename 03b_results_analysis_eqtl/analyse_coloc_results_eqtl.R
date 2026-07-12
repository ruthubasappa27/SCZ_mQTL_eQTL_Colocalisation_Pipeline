# SCZ mQTL eQTL Colocalisation Pipeline
# Script: eQTL Results Analysis and Per-Gene Colocalisation
# Author: Ruthushree P Basappa
# Institution: University of Galway
# Supervisor: Prof. Derek Morris
# Date: June 2026

library(data.table)
library(dplyr)
library(here)

here::i_am("03b_results_analysis_eqtl/analyse_coloc_results_eqtl.R")

dir.create(here("results_eqtl"), showWarnings = FALSE, recursive = TRUE)

# Step 1: Load locus-level colocalisation results

all_results <- rbindlist(
  lapply(
    list.files(
      here("results_eqtl"),
      pattern = "results_summary_coloc_abf.csv",
      recursive = TRUE,
      full.names = TRUE
    ),
    fread
  ),
  fill = TRUE
)

all_results_unique <- all_results[!duplicated(region)]
cat("Total unique loci tested:", nrow(all_results_unique), "\n")
cat("Loci passing PP4 > 0.8 (locus level):", sum(all_results_unique$PP.H4.abf > 0.8, na.rm = TRUE), "\n")

locus_hits <- all_results_unique[PP.H4.abf > 0.8][order(-PP.H4.abf)]
print(locus_hits[, .(region, nsnps, PP.H4.abf)])

fwrite(locus_hits, here("results_eqtl", "coloc_eqtl_hits.csv"))

# Step 2: Per-gene colocalisation at loci passing the locus-level threshold
# Following Dobbyn et al. (2018): for each locus passing PP4 >= 0.8, coloc.abf
# was run independently for each gene with available eQTL data within that
# locus window, to identify which specific gene(s) drove the colocalised signal.

pergene_results <- rbindlist(
  lapply(
    list.files(
      here("results_eqtl"),
      pattern = "coloc_pergene_.*\\.csv",
      recursive = TRUE,
      full.names = TRUE
    ),
    fread
  ),
  fill = TRUE
)

fwrite(pergene_results, here("results_eqtl", "coloc_eqtl_pergene_results.csv"))

pergene_significant <- pergene_results[PP.H4.abf >= 0.8][order(-PP.H4.abf)]
cat("Gene-level signals passing PP4 >= 0.8:", nrow(pergene_significant), "\n")
print(pergene_significant[, .(hgnc_symbol, region, nsnps, PP.H4.abf)])

fwrite(pergene_significant, here("results_eqtl", "coloc_eqtl_pergene_significant.csv"))

# Step 3: Directionality analysis for significant gene-level signals

directionality <- data.table()
for (i in seq_len(nrow(pergene_significant))) {
  gene <- pergene_significant$hgnc_symbol[i]
  region <- pergene_significant$region[i]

  harm_file <- here("results_eqtl", paste0(gsub("[:-]", "_", region), "_coloc"),
                     "data", "datasets", "harmonised_sumstats.csv")

  if (!file.exists(harm_file)) next

  harm <- fread(harm_file)
  top_snp <- harm[trait == "SCZ"][order(pvalues)][1, snp]

  dir_result <- harm %>%
    dplyr::filter(snp == top_snp) %>%
    dplyr::select(snp, trait, beta, pvalues)

  scz_beta <- dir_result$beta[dir_result$trait == "SCZ"]
  eqtl_beta <- dir_result$beta[dir_result$trait == "eQTL"]

  if (length(scz_beta) > 0 && length(eqtl_beta) > 0) {
    direction <- ifelse(sign(scz_beta) == sign(eqtl_beta), "Increases", "Decreases")
    directionality <- rbind(directionality,
                             data.table(Gene = gene, TopSNP = top_snp,
                                        GWAS_beta = scz_beta, eQTL_beta = eqtl_beta,
                                        Direction = direction))
  }
}

fwrite(directionality, here("results_eqtl", "coloc_eqtl_directionality.csv"))
cat("Directionality analysis saved for", nrow(directionality), "genes\n")

# Step 4: Credible set construction for significant gene-level signals

credible_sets <- data.table()
for (i in seq_len(nrow(pergene_significant))) {
  gene <- pergene_significant$hgnc_symbol[i]
  region <- pergene_significant$region[i]

  snpwise_file <- here("results_eqtl", paste0(gsub("[:-]", "_", region), "_coloc"),
                        "tables", "coloc.abf_snpwise_PP_H4_abf.csv")

  if (!file.exists(snpwise_file)) next

  snpwise <- fread(snpwise_file)
  snpwise <- snpwise[order(-SNP.PP.H4)]
  snpwise[, cumPP := cumsum(SNP.PP.H4) / sum(SNP.PP.H4)]
  cs <- snpwise[cumPP <= 0.95]
  cs[, Gene := gene]

  credible_sets <- rbind(credible_sets, cs, fill = TRUE)
}

fwrite(credible_sets, here("results_eqtl", "coloc_eqtl_credible_sets_full.csv"))
cat("Credible sets saved for", uniqueN(credible_sets$Gene), "genes\n")
