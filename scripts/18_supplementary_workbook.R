# 18_supplementary_workbook.R
# Supplementary Tables S1-S9 in one Excel workbook, numbered as in the manuscript.
# Values are written unrounded; Excel number formats control how they are displayed.
# Missing values: NA (not applicable or not computed); "none" where a CpG site has no gene annotation.
# Output: results/supplementary/Supplementary_Tables.xlsx
# Usage:  Rscript scripts/18_supplementary_workbook.R
# Requires openxlsx.

source("config/config.R")
source("scripts/functions.R")
if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package openxlsx is required: install.packages(\"openxlsx\")")
start_log("18_supplementary_workbook")
res <- PATHS$results

read_tsv <- function(...) {
  f <- file.path(res, ...)
  if (!file.exists(f)) stop("Missing ", f)
  fread(f)
}

# "A;A;B" -> "A;B"; empty -> fill
tidy_genes <- function(x, fill) {
  x <- vapply(strsplit(ifelse(is.na(x), "", x), ";"), function(g) paste(unique(trimws(g[g != ""])), collapse = ";"), "")
  ifelse(x == "", fill, x)
}

# spec: one row per output column - source column, header, number format
# formats: text, int, pos, pp (3 decimals), num (4 decimals), mean (2 decimals), p (scientific), yesno
build <- function(dt, spec) {
  spec <- as.data.table(do.call(rbind, spec))
  setnames(spec, c("src", "header", "fmt"))
  missing <- setdiff(spec$src, names(dt))
  if (length(missing)) stop("Columns not found: ", paste(missing, collapse = ", "))
  out <- dt[, spec$src, with = FALSE]
  for (i in which(spec$fmt == "yesno")) {
    v <- out[[i]]
    set(out, j = i, value = ifelse(is.na(v), NA_character_, ifelse(as.logical(v), "Yes", "No")))
  }
  for (i in which(spec$fmt == "text")) {
    v <- as.character(out[[i]])
    set(out, j = i, value = ifelse(v == "", NA_character_, v))
  }
  setnames(out, spec$header)
  list(data = out, fmt = spec$fmt)
}

# S1: quality control of the input datasets and locus definition (steps 01-04)
qc <- read_tsv("qc_counts.tsv")
qc_labels <- c("01_qc_gwas" = "GWAS", "03_qc_mqtl" = "mQTL", "04_qc_eqtl" = "eQTL", "02_define_loci" = "Loci")
if (length(setdiff(names(qc_labels), qc$script)))
  stop("No counts in qc_counts.tsv for: ", paste(setdiff(names(qc_labels), qc$script), collapse = ", "))
qc <- qc[script %in% names(qc_labels)][, dataset := qc_labels[script]]
s1 <- build(qc, list(c("dataset", "Dataset", "text"), c("step", "Step", "text"), c("n", "Count", "int")))

# S2: loci
loci <- read_tsv("loci", "loci_all.tsv")
loci[, status := fifelse(status == "analysed", "Analysed", fifelse(status == "excluded_MHC", "Excluded (MHC)", status))]
s2 <- build(loci, list(
  c("locus_id", "Locus", "text"), c("CHR", "Chromosome", "int"),
  c("start", "Start (GRCh37)", "pos"), c("end", "End (GRCh37)", "pos"),
  c("lead_snp", "Lead SNP", "text"), c("lead_p", "Lead SNP p", "p"),
  c("n_lead_snps", "Lead SNPs merged", "int"), c("merged_leads", "Merged lead SNPs", "text"),
  c("status", "Status", "text")))

# S3/S4: colocalisation tests, with details of colocalised features
coloc_data <- function(qtl) {
  keys <- c("locus_id", "feature")
  out <- read_tsv(paste0("coloc_", qtl), "coloc_all_tests.tsv")[status == "tested"]
  out[, row_order := .I]
  for (f in c("hits.tsv", "hits_directionality.tsv", "hits_gene_annotation.tsv")) {
    add <- read_tsv(paste0("coloc_", qtl), f)
    extra <- setdiff(names(add), names(out))
    if (length(extra)) out <- merge(out, add[, c(keys, extra), with = FALSE], by = keys, all.x = TRUE, sort = FALSE)
  }
  out[, colocalised := !is.na(effect_of_risk_allele)]
  setorder(out, row_order)[, row_order := NULL]
  out
}
coloc_common <- list(
  c("nsnps", "SNPs tested", "int"), c("n_freq_mismatch", "SNPs removed (allele frequency)", "int"),
  c("gwas_minp_overlap", "Minimum GWAS p", "p"), c("qtl_minp_overlap", "Minimum QTL p", "p"),
  c("gwas_lead_in_overlap", "GWAS lead SNP tested", "yesno"),
  c("method", "Primary method", "text"), c("reason", "Reason for primary method", "text"),
  c("PP.H0", "PP0", "pp"), c("PP.H1", "PP1", "pp"), c("PP.H2", "PP2", "pp"), c("PP.H3", "PP3", "pp"),
  c("PP4", "PP4 (primary method)", "pp"),
  c("abf_PP.H4", "PP4 coloc.abf", "pp"), c("susie_PP.H4_max", "PP4 coloc.susie", "pp"),
  c("susie_ncs_gwas", "SuSiE credible sets, GWAS", "int"), c("susie_ncs_qtl", "SuSiE credible sets, QTL", "int"),
  c("susie_pair", "SuSiE credible-set pair", "text"),
  c("PP4_p12_1e-06", "PP4, p12 = 1e-6", "pp"), c("PP4_p12_5e-06", "PP4, p12 = 5e-6", "pp"),
  c("PP4_p12_1e-05", "PP4, p12 = 1e-5", "pp"), c("PP4_p12_5e-05", "PP4, p12 = 5e-5", "pp"),
  c("top_snp", "Top SNP", "text"), c("top_snp_pp", "Top SNP posterior probability", "pp"),
  c("cs95_size", "95% credible set size", "int"),
  c("ld_mismatch_s_gwas", "LD mismatch s, GWAS", "num"), c("ld_mismatch_s_qtl", "LD mismatch s, QTL", "num"),
  c("robust_to_p12", "Stable across priors", "yesno"),
  c("gwas_lead", "GWAS lead SNP", "text"), c("r2_top_vs_gwas_lead", "r2, top SNP and GWAS lead", "num"),
  c("cs_max_r2_gwas_lead", "Maximum r2, credible set and GWAS lead", "num"),
  c("risk_allele", "Risk allele", "text"), c("other_allele", "Other allele", "text"),
  c("gwas_OR_risk", "GWAS odds ratio (risk allele)", "num"), c("gwas_p", "GWAS p (top SNP)", "p"),
  c("qtl_beta_risk", "QTL effect (risk allele)", "num"), c("qtl_p", "QTL p (top SNP)", "p"),
  c("effect_of_risk_allele", "Direction", "text"))

m <- coloc_data("mqtl")
m[, feature_label := tidy_genes(feature_label, "none")]
m[colocalised == TRUE, gene_symbol := tidy_genes(gene_symbol, "none")]
s3 <- build(m, c(list(
  c("locus_id", "Locus", "text"), c("feature", "CpG site", "text"),
  c("feature_chr", "Chromosome", "int"), c("feature_pos", "Position (GRCh37)", "pos"),
  c("feature_label", "Gene (Illumina)", "text")), coloc_common, list(
  c("gene_symbol", "Gene (HGNC)", "text"),
  c("nearest_gene", "Nearest gene", "text"), c("nearest_gene_dist_bp", "Distance to nearest gene (bp)", "pos"),
  c("nearest_coding", "Nearest protein-coding gene", "text"),
  c("nearest_coding_dist_bp", "Distance to nearest protein-coding gene (bp)", "pos"),
  c("genes_within_50kb", "Genes within 50 kb", "text"),
  c("illumina_agrees_nearest", "Illumina gene is nearest gene", "yesno"))))

e <- coloc_data("eqtl")
s4 <- build(e, c(list(
  c("locus_id", "Locus", "text"), c("feature", "Gene (Ensembl)", "text"),
  c("feature_chr", "Chromosome", "int"), c("feature_pos", "TSS (GRCh38)", "pos")),
  coloc_common, list(c("gene_symbol", "Gene symbol", "text"), c("biotype", "Biotype", "text"))))

coloc_note <- paste("NA, not applicable: SuSiE columns where coloc.susie was not run; columns from",
                    "'Stable across priors' onwards for features that did not colocalise (PP4 \u2265 0.8).")

# S5/S6: SMR and HEIDI
smr_note_text <- c("probe not in BESD" = "Not in SMR-format data",
                   "no instrument at this p threshold" = "No instrument at this threshold",
                   "HEIDI not computable (inconclusive)" = "HEIDI not computable")
smr_common <- list(
  c("analysis", "Analysis", "text"), c("peqtl", "Instrument threshold (p)", "p"),
  c("coloc_hits", "Colocalised locus (PP4)", "text"),
  c("topSNP", "Top SNP", "text"), c("topSNP_chr", "Chromosome", "int"), c("topSNP_bp", "Top SNP position", "pos"),
  c("A1", "Effect allele", "text"), c("A2", "Other allele", "text"), c("Freq", "Effect allele frequency", "pp"),
  c("b_GWAS", "b GWAS", "num"), c("se_GWAS", "SE GWAS", "num"), c("p_GWAS", "p GWAS", "p"),
  c("b_eQTL", "b QTL", "num"), c("se_eQTL", "SE QTL", "num"), c("p_eQTL", "p QTL", "p"),
  c("b_SMR", "b SMR", "num"), c("se_SMR", "SE SMR", "num"), c("p_SMR", "p SMR", "p"),
  c("n_tested", "Tests (Bonferroni)", "int"), c("smr_bonferroni_sig", "Significant after Bonferroni", "yesno"),
  c("p_HEIDI", "p HEIDI", "p"), c("nsnp_HEIDI", "SNPs in HEIDI", "int"),
  c("heidi_pass", "Passed HEIDI", "yesno"), c("note", "Note", "text"))
smr_prep <- function(dt) {
  dt[, analysis := fifelse(analysis == "primary", "Primary", fifelse(analysis == "sensitivity", "Sensitivity", analysis))]
  dt[note %in% names(smr_note_text), note := smr_note_text[note]]
  dt
}

sm <- smr_prep(read_tsv("smr_mqtl", "smr_results.tsv"))
sm[, feature_label := tidy_genes(feature_label, "none")]
s5 <- build(sm, c(list(c("probeID", "CpG site", "text"), c("feature_label", "Gene (Illumina)", "text")), smr_common))

se <- smr_prep(read_tsv("smr_eqtl", "smr_results.tsv"))[, row_order := .I]
sym <- read_tsv("coloc_eqtl", "hits_gene_annotation.tsv")[, .(probeID = feature, gene_symbol)]
se <- setorder(merge(se, unique(sym), by = "probeID", all.x = TRUE, sort = FALSE), row_order)
s6 <- build(se, c(list(c("probeID", "Gene (Ensembl)", "text"), c("gene_symbol", "Gene symbol", "text")), smr_common))
smr_note <- "NA, not available; the reason is given in the Note column."

# S7: developmental regulation
s7 <- build(read_tsv("characterisation", "dev_regulation.tsv"), list(
  c("row", "Gene", "text"),
  c("n_fetal_brains", "Fetal brains", "int"), c("n_postnatal_brains", "Postnatal brains", "int"),
  c("fetal_mean", "Mean log2 expression, fetal", "mean"), c("postnatal_mean", "Mean log2 expression, postnatal", "mean"),
  c("log2_diff_fetal_vs_postnatal", "log2 difference, fetal vs postnatal", "mean"),
  c("p", "p (Welch t-test)", "p"), c("fdr", "FDR", "p"),
  c("max_mean_log2", "Maximum mean log2 expression", "mean"),
  c("expressed", "Expressed", "yesno"), c("direction", "Higher in", "text"),
  c("regulated", "Developmentally regulated", "yesno")))

# S8: cell-type expression
s8 <- build(read_tsv("characterisation", "celltype_summary.tsv"), list(
  c("current_symbol", "Gene", "text"), c("in_codex", "In CoDEx", "yesno"),
  c("peak_cell_type", "Cell type with highest expression", "text"),
  c("peak_mean", "Highest mean expression", "pp"), c("min_mean", "Lowest mean expression", "pp"),
  c("range", "Range", "pp")))

# S9: credible sets
cs <- rbindlist(lapply(c(mQTL = "mqtl", eQTL = "eqtl"), function(q)
  read_tsv(paste0("coloc_", q), "hits_credible_sets.tsv")), idcol = "dataset")
cs[dataset == "mQTL", feature_label := tidy_genes(feature_label, "none")]
s9 <- build(cs, list(
  c("dataset", "Dataset", "text"), c("locus_id", "Locus", "text"),
  c("feature", "CpG site or gene", "text"), c("feature_label", "Gene annotation", "text"),
  c("snp", "SNP", "text"), c("pp", "Posterior probability", "pp"), c("cum_pp", "Cumulative posterior probability", "pp"),
  c("r2_with_top_snp", "r2 with top SNP", "num"), c("r2_with_gwas_lead", "r2 with GWAS lead", "num"),
  c("nearest_gwas_lead", "Nearest GWAS lead SNP", "text")))

tables <- list(
  list("S1 QC", "Quality-control counts for each dataset.", s1),
  list("S2 Loci", "Schizophrenia GWAS loci.", s2),
  list("S3 mQTL colocalisation", "All mQTL colocalisation tests, with PP4 from both methods, prior sensitivity, SNPs removed by the allele-frequency check, and direction of effect and gene annotation of colocalised CpG sites.",
       s3, paste(coloc_note, "none, no gene annotation.")),
  list("S4 eQTL colocalisation", "All eQTL colocalisation tests, with PP4 from both methods, prior sensitivity, SNPs removed by the allele-frequency check, and direction of effect and gene annotation of colocalised genes.",
       s4, coloc_note),
  list("S5 SMR HEIDI CpG sites", "SMR and HEIDI results for the colocalised CpG sites at the primary and sensitivity instrument thresholds.",
       s5, paste(smr_note, "none, no gene annotation.")),
  list("S6 SMR HEIDI genes", "SMR and HEIDI results for the colocalised genes at the primary and sensitivity instrument thresholds.",
       s6, smr_note),
  list("S7 Developmental expression", "Developmental regulation of genes annotated to the prioritised CpG sites.", s7),
  list("S8 Cell-type expression", "Cell-type expression of genes annotated to the prioritised CpG sites.",
       s8, "NA, gene not present in CoDEx."),
  list("S9 Credible sets", "95% credible sets of the colocalised CpG sites and genes.",
       s9, "none, no gene annotation.")
)

# workbook
wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 10, fontName = "Calibri")
bold <- openxlsx::createStyle(textDecoration = "bold")
header <- openxlsx::createStyle(textDecoration = "bold", wrapText = TRUE, valign = "top", border = "bottom")
num_fmt <- c(int = "0", pos = "0", pp = "0.000", num = "0.0000", mean = "0.00", p = "0.00E+00")
col_width <- 16

openxlsx::addWorksheet(wb, "Contents")
contents <- data.table(Table = paste0("Supplementary Table ", sub(" .*", "", sapply(tables, `[[`, 1))),
                       Sheet = sapply(tables, `[[`, 1),
                       Description = sapply(tables, `[[`, 2),
                       Rows = sapply(tables, function(t) nrow(t[[3]]$data)))
openxlsx::writeData(wb, "Contents", contents, headerStyle = bold)
openxlsx::setColWidths(wb, "Contents", cols = 1:4, widths = c(24, 28, 110, 8))

for (t in tables) {
  sheet <- t[[1]]; d <- t[[3]]$data; fmt <- t[[3]]$fmt
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, paste0("Supplementary Table ", sub(" .*", "", sheet), ". ", t[[2]]),
                      startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet, bold, rows = 1, cols = 1)
  if (length(t) >= 4) openxlsx::writeData(wb, sheet, t[[4]], startRow = 2, startCol = 1)
  openxlsx::writeData(wb, sheet, d, startRow = 3, headerStyle = header, keepNA = TRUE, na.string = "NA")
  rows <- 3 + seq_len(nrow(d))
  for (f in names(num_fmt)) {
    cols <- which(fmt == f)
    if (length(cols) && nrow(d))
      openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = num_fmt[[f]]),
                         rows = rows, cols = cols, gridExpand = TRUE)
  }
  openxlsx::freezePane(wb, sheet, firstActiveRow = 4)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(d)), widths = col_width)
  openxlsx::setRowHeights(wb, sheet, rows = 3, heights = 40)
  logmsg(sprintf("  %-60s %s", sheet, format(nrow(d), big.mark = ",")))
}

out <- file.path(res, "supplementary", "Supplementary_Tables.xlsx")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
logmsg("  Written: ", out)
finish_log()
