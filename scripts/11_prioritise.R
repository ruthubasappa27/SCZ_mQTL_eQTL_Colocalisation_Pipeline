# 11_prioritise.R
# Prioritised features: colocalised (PP4 >= 0.8), SMR significant after Bonferroni correction at
# the primary instrument threshold, and HEIDI p > 0.05. Gene names are the current HGNC symbols
# from 08_annotate_genes.R; multi-gene CpG labels are split at ";" and unlabelled CpG sites
# contribute no gene. eQTL genes without an HGNC symbol keep their Ensembl ID.
#
# Output: results/characterisation/prioritised_genes.tsv
# Usage:  Rscript scripts/11_prioritise.R

source("config/config.R")
source("scripts/functions.R")
start_log("11_prioritise")
outdir <- file.path(PATHS$results, "characterisation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

prior <- function(qtl) {
  h <- fread(file.path(PATHS$results, paste0("coloc_", qtl), "hits.tsv"))
  s <- fread(file.path(PATHS$results, paste0("smr_", qtl), "smr_results.tsv"))[analysis == "primary"]
  a <- fread(file.path(PATHS$results, paste0("coloc_", qtl), "hits_gene_annotation.tsv"))
  s <- s[smr_bonferroni_sig == TRUE & heidi_pass %in% TRUE]
  h <- h[feature %in% s$probeID]
  if (qtl == "mqtl") {
    lab <- setNames(as.character(a$gene_symbol), a$feature)
    g <- h[, .(gene = trimws(unlist(strsplit(ifelse(is.na(lab[feature]), "", lab[feature]), ";")))),
           by = .(feature, locus_id)]
  } else {
    sym <- setNames(as.character(a$gene_symbol), a$feature)
    g <- h[, .(gene = ifelse(is.na(sym[feature]) | sym[feature] == "", feature, sym[feature])),
           by = .(feature, locus_id)]
  }
  count_step(sprintf("%s prioritised features", qtl), nrow(h))
  count_step(sprintf("%s prioritised features without a gene name", qtl),
             nrow(h) - uniqueN(g[gene != ""]$feature))
  g[gene != "", .(qtl, feature, locus_id, gene)]
}
genes_dt <- rbind(prior("mqtl"), prior("eqtl"))
fwrite(genes_dt, file.path(outdir, "prioritised_genes.tsv"), sep = "\t")
count_step("Unique gene names (mQTL + eQTL combined)", uniqueN(genes_dt$gene))
finish_log()
