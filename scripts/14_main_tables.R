# 14_main_tables.R
# Table 2: prioritised CpG sites (colocalised PP4 >= 0.8, SMR significant at the primary
# threshold, HEIDI p > 0.05), one row per CpG site, ordered by genomic position. Gene names are
# the current HGNC symbols from 08_annotate_genes.R. PP4 is given for coloc.susie and coloc.abf,
# to six decimal places; the value from the primary method is in bold.
#
# Output: results/tables/Table2.docx, Table2.tsv
# Usage:  Rscript scripts/14_main_tables.R
# Requires flextable and officer.

source("config/config.R")
source("scripts/functions.R")
suppressMessages({ library(flextable); library(officer) })
start_log("14_main_tables")
outdir <- file.path(PATHS$results, "tables")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fmt_pp <- function(x) ifelse(is.na(x), "n.e.", formatC(x, format = "f", digits = 6))
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
compress <- function(g) {   # PCDHA1, PCDHA2, PCDHA3, PCDHA4 -> PCDHA1-4
  fam <- sub("[0-9]+$", "", g); num <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", g, perl = TRUE)))
  out <- character()
  for (f in unique(fam)) {
    k <- which(fam == f)
    if (length(k) > 2 && all(!is.na(num[k])) && all(diff(sort(num[k])) == 1))
      out <- c(out, sprintf("%s%d-%d", f, min(num[k]), max(num[k]))) else out <- c(out, g[k])
  }
  out
}

build <- function(q) {
  h <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits.tsv"))
  s <- fread(file.path(PATHS$results, paste0("smr_", q), "smr_results.tsv"))[analysis == "primary"]
  a <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits_gene_annotation.tsv"))
  dr <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits_directionality.tsv"))
  g <- fread(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"), select = c("SNP", "CHR", "POS"))
  s <- s[smr_bonferroni_sig == TRUE & heidi_pass %in% TRUE]
  d <- h[feature %in% s$probeID, .(locus_id, feature, PP4, method, top_snp,
                                    pp4_susie = fifelse(method == "coloc.susie", PP4, susie_PP.H4_max),
                                    pp4_abf = abf_PP.H4)]
  d <- merge(d, s[, .(feature = probeID, p_SMR, p_HEIDI)], by = "feature")
  d <- merge(d, dr[, .(feature, risk_allele, effect_of_risk_allele)], by = "feature")
  d <- merge(d, a[, .(feature, gene_symbol)], by = "feature", all.x = TRUE)
  d <- merge(d, g, by.x = "top_snp", by.y = "SNP")
  if (nrow(d) != nrow(s)) stop(q, ": rows lost when combining result files")
  d[, gene := vapply(strsplit(fifelse(is.na(gene_symbol), "", gene_symbol), ";"), function(x) {
    x <- x[nzchar(x) & !grepl("^ENSG", x)]; if (length(x)) paste(compress(x), collapse = ", ") else "none" }, "")]
  setorder(d, CHR, POS, feature)
  d[, `:=`(snp_allele = sprintf("%s (%s)", top_snp, risk_allele),
           susie_txt  = fmt_pp(pp4_susie),
           abf_txt    = fmt_pp(pp4_abf),
           smr_heidi  = sprintf("%s / %s", fmt_p(p_SMR), fmt_p(p_HEIDI)),
           effect     = fifelse(grepl("^increases", effect_of_risk_allele), "\u2191", "\u2193"))]
  count_step(sprintf("%s table rows", q), nrow(d))
  d
}

m <- build("mqtl")
fwrite(m[, .(Gene = gene, CpG = feature, `Top SNP (risk allele)` = snp_allele, `Primary method` = method,
             `PP4 SuSiE` = pp4_susie, `PP4 ABF` = pp4_abf, p_SMR, p_HEIDI, Effect = effect)],
       file.path(outdir, "Table2.tsv"), sep = "\t")

style <- function(ft) {
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = fp_border(width = 1), part = "header")
  ft <- hline_bottom(ft, border = fp_border(width = 0.75), part = "header")
  ft <- hline_bottom(ft, border = fp_border(width = 1), part = "body")
  ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
  autofit(ft)
}
bold_primary <- function(ft, d) {
  for (x in list(c("coloc.susie", "PP4 SuSiE"), c("coloc.abf", "PP4 ABF"))) {
    i <- which(d$method == x[1])
    if (length(i)) ft <- bold(ft, i = i, j = x[2], part = "body")
  }
  ft
}

# Table 2: mQTL; the gene name is shown once for consecutive CpGs of the same gene
t2 <- m[, .(Gene = gene, CpG = feature, `Top SNP (risk allele)` = snp_allele,
            `PP4 SuSiE` = susie_txt, `PP4 ABF` = abf_txt, `SMR p / HEIDI p` = smr_heidi, Effect = effect)]
ft2 <- if (nrow(m)) valign(bold_primary(style(merge_v(flextable(t2), j = "Gene")), m), j = "Gene", valign = "top", part = "body")

note <- paste0("PP4, posterior probability of a shared causal variant from coloc.susie (SuSiE) and coloc.abf (ABF); ",
               "bold, primary method (SuSiE when both traits have a credible set, otherwise ABF); ",
               "n.e., not estimable (no GWAS credible set). ",
               "SMR p, Bonferroni-significant in all rows; HEIDI p > 0.05 in all rows. ",
               "\u2191/\u2193, effect of the schizophrenia risk allele on ")
add_table <- function(doc, title, ft, n, footnote, what) {
  doc <- body_add_par(doc, title, style = "Normal")
  if (n == 0) return(body_add_par(doc, sprintf("No %s passed colocalisation, SMR and HEIDI.", what), style = "Normal"))
  doc <- body_add_flextable(doc, ft)
  body_add_par(doc, footnote, style = "Normal")
}
doc <- read_docx()
doc <- add_table(doc, "Table 2. CpG sites prioritised by colocalisation, SMR and HEIDI", ft2, nrow(m),
                 paste0(note, "DNA methylation. none, no gene annotation."), "mQTL CpG")
print(doc, target = file.path(outdir, "Table2.docx"))
logmsg("  Table written to ", file.path(outdir, "Table2.docx"))
finish_log()
