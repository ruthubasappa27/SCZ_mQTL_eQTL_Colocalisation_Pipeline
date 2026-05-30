library(data.table)
library(dplyr)

# Reload everything
gwas <- fread("~/Desktop/Thesis/gwas_qc.csv")
mqtl_full <- fread("~/Desktop/Thesis/mqtl_preprocessed.csv")

all_results <- rbindlist(lapply(
  list.files("~/Desktop/Thesis/COLOC-reporter/coloc/results", 
             pattern="results_summary_coloc_abf.csv", 
             recursive=TRUE, full.names=TRUE),
  fread), fill=TRUE)

all_results_unique <- all_results[!duplicated(region)]
hits <- all_results_unique[PP.H4.abf > 0.8][order(-PP.H4.abf)]