# 12_developmental_expression.R
# Developmental expression of the prioritised genes in the Kang et al. (2011) human brain
# transcriptome (GEO GSE25219, platform GPL5175, transcript-cluster level, log2 signal).
#   Genes:   results/characterisation/prioritised_genes.tsv (11_prioritise.R), matched
#            to GPL5175 transcript clusters by Entrez Gene ID (org.Hs.eg.db; eQTL genes from the
#            Ensembl ID, mQTL genes from symbol or alias). Genes are labelled with their current
#            symbol; transcript clusters of one gene are averaged.
#   Regions: NCX (11 neocortical areas, and frontal, parietal, temporal, occipital and
#            M1C&S1C cortex in early fetal samples), HIP, AMY, STR, MD, CBC. Transient fetal
#            structures (e.g. ganglionic eminences, upper rhombic lip) are listed and not used.
#   Periods: 1-15, from the "Stage" field (1-7 prenatal, 8-15 postnatal).
#   Developmental regulation (thresholds of Kang et al. 2011): expressed = mean log2 >= 6 in at
#            least one region and period; regulated = fetal (periods 3-7) vs postnatal (8-12)
#            difference >= 1 log2 (2-fold) at BH FDR < 0.01, in either direction; each brain is
#            the mean of its regions (one value per brain), compared with a Welch t-test; FDR
#            is computed separately for the mQTL and eQTL genes.
#   Heatmaps: mean log2 expression, all periods; cells with fewer than 3 samples are grey; genes
#            measured by the same transcript clusters share one row. Trajectory plots of the
#            developmentally regulated rows. mQTL and eQTL genes are plotted separately; plot
#            objects are saved (.rds) for 17_main_figures.R.
#
# Output: results/characterisation/dev_expression_*.tsv, dev_regulation.tsv,
#         dev_expression_heatmap_<qtl>.pdf/.png/.rds, dev_trajectories_regulated_<qtl>.pdf/.png/.rds
# Usage:  Rscript scripts/12_developmental_expression.R
# Figures are drawn at printed size (17 cm wide; text >= 7 pt).
# Requires GEOquery and org.Hs.eg.db (Bioconductor) and ggplot2; downloads GSE25219 on first use.

source("config/config.R")
source("scripts/functions.R")
suppressMessages({ library(GEOquery); library(ggplot2); library(org.Hs.eg.db) })
start_log("12_developmental_expression")
outdir <- file.path(PATHS$results, "characterisation")
refdir <- PATHS$kang_dir
dir.create(refdir, showWarnings = FALSE, recursive = TRUE)

ncx_areas <- c("OFC", "DFC", "VFC", "MFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C",
               "FC", "PC", "TC", "OC", "M1C&S1C")   # early fetal cortex sampled by lobe
other_regions <- c("HIP", "AMY", "STR", "MD", "CBC")

# genes
pg <- fread(file.path(outdir, "prioritised_genes.tsv"))
genes <- unique(pg[, .(gene, qtl, ensembl = fifelse(qtl == "eqtl", sub("\\..*$", "", feature), NA_character_))])
count_step("Prioritised genes", nrow(genes))

# data
mfile <- file.path(refdir, "GSE25219-GPL5175_series_matrix.txt.gz")
if (!file.exists(mfile)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25219/matrix/GSE25219-GPL5175_series_matrix.txt.gz"
  logmsg("  Downloading ", url, " (large; one time only)")
  if (system2("curl", c("-sSfL", "-o", shQuote(paste0(mfile, ".part")), shQuote(url))) != 0)
    stop("Download failed: ", url, ".")
  invisible(file.rename(paste0(mfile, ".part"), mfile))
}
es  <- getGEO(filename = mfile, getGPL = FALSE)
ex  <- exprs(es)
ph  <- as.data.table(pData(es), keep.rownames = "sample")
count_step("GSE25219 samples", ncol(ex))
count_step("GSE25219 transcript clusters", nrow(ex))
if (max(ex, na.rm = TRUE) > 30) stop("Values do not look like log2 signal (max ", max(ex, na.rm = TRUE), ").")

# region and period columns in the sample annotation
char_col <- function(key) {
  cols <- grep("characteristics", names(ph), value = TRUE)
  hit <- cols[vapply(cols, function(c) any(grepl(paste0("^", key, ":"), ph[[c]], ignore.case = TRUE)), logical(1))]
  if (length(hit) == 0) {
    print(head(ph[, ..cols], 3))
    stop("No '", key, ":' field in the GSE25219 sample annotation (fields printed above).")
  }
  trimws(sub("^[^:]*:", "", ph[[hit[1]]]))
}
ph[, region := char_col("region")]
ph[, period := as.integer(char_col("stage"))]
ph[, donor := char_col("brain code")]
logmsg("  Samples per region code:")
print(ph[, .N, by = region][order(region)])
ph[, region6 := fifelse(region %in% ncx_areas, "NCX", fifelse(region %in% other_regions, region, NA_character_))]
count_step("Samples in NCX, HIP, AMY, STR, MD or CBC", sum(!is.na(ph$region6)))
count_step("Samples with another region code (not used)", sum(is.na(ph$region6)))

gpl <- getGEO("GPL5175", destdir = refdir)
ann <- as.data.table(Table(gpl))[, .(ID = as.character(ID), gene_assignment)]
ann <- ann[!is.na(gene_assignment) & gene_assignment != "---"]
ann <- ann[, {
  p <- strsplit(unlist(strsplit(gene_assignment, " /// ")), " // ")
  .(symbol = vapply(p, function(x) if (length(x) >= 2) trimws(x[2]) else NA_character_, ""),
    entrez = vapply(p, function(x) if (length(x) >= 5) trimws(x[5]) else NA_character_, ""))
}, by = ID]
ann <- unique(ann[!is.na(entrez) & entrez != "---"])

# gene -> Entrez Gene ID -> transcript clusters
# eQTL genes by Ensembl ID; mQTL genes by symbol, then by alias. Only unique mappings are used.
eg_of <- function(keys, keytype) {
  if (!length(keys)) return(character())
  r <- tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db, keys, "ENTREZID", keytype, multiVals = "list"),
                error = function(e) setNames(vector("list", length(keys)), keys))
  vapply(keys, function(k) { x <- unique(na.omit(unlist(r[[k]]))); if (length(x) == 1) x else NA_character_ }, "")
}
genes[, entrez := NA_character_][, match := NA_character_]
genes[!is.na(ensembl), entrez := eg_of(ensembl, "ENSEMBL")]
genes[!is.na(entrez) & is.na(match), match := "Ensembl ID"]
genes[is.na(entrez), entrez := eg_of(gene, "SYMBOL")]
genes[!is.na(entrez) & is.na(match), match := "symbol"]
genes[is.na(entrez), entrez := eg_of(gene, "ALIAS")]
genes[!is.na(entrez) & is.na(match), match := "alias"]
genes[, current_symbol := fifelse(is.na(entrez), gene,
  unname(AnnotationDbi::mapIds(org.Hs.eg.db, fifelse(is.na(entrez), "1", entrez), "SYMBOL", "ENTREZID")))]
genes[, name_used := vapply(entrez, function(e) if (is.na(e)) NA_character_ else
  paste(unique(ann[entrez == e]$symbol), collapse = "/"), "")]
genes[name_used == "", name_used := NA_character_]
genes[, n_transcript_clusters := vapply(entrez, function(e) if (is.na(e)) 0L else
  length(intersect(ann[entrez == e]$ID, rownames(ex))), 0L)]
genes[, found := n_transcript_clusters > 0]
count_step("Genes with a unique Entrez Gene ID", sum(!is.na(genes$entrez)))
count_step("Genes found in GSE25219", sum(genes$found))
count_step("Genes not found (reported as not available)", sum(!genes$found))

# mean log2 expression per gene, region and period
keep <- ph[!is.na(region6) & !is.na(period)]
long <- rbindlist(lapply(which(genes$found), function(i) {
  ids <- intersect(ann[entrez == genes$entrez[i]]$ID, rownames(ex))
  v <- colMeans(ex[ids, keep$sample, drop = FALSE], na.rm = TRUE)
  keep[, .(gene = genes$current_symbol[i], name_used = genes$name_used[i], qtl = genes$qtl[i],
           n_transcript_clusters = length(ids), region = region6, period, value = v)][
    , .(mean_log2 = mean(value, na.rm = TRUE), n_samples = .N),
    by = .(gene, name_used, qtl, n_transcript_clusters, region, period)]
}))
fwrite(genes, file.path(outdir, "dev_expression_gene_availability.tsv"), sep = "\t")
fwrite(long[order(qtl, gene, region, period)], file.path(outdir, "dev_expression_values.tsv"), sep = "\t")

# rows: genes measured by the same transcript clusters are combined
genes[found == TRUE, tc_key := vapply(entrez, function(e)
  paste(sort(intersect(ann[entrez == e]$ID, rownames(ex))), collapse = ","), "")]
rows <- genes[found == TRUE, .(row = paste(sort(unique(current_symbol)), collapse = "/"),
                               arm = paste(sort(unique(qtl)), collapse = "/"),
                               entrez = entrez[1]), by = tc_key]
gm <- unique(genes[found == TRUE, .(gene = current_symbol, tc_key)][rows, on = "tc_key"][, .(gene, row, arm)])
dev <- long[, .(gene, region, period, mean_log2, n_samples)]
dev <- unique(merge(dev, gm, by = "gene")[, .(row, arm, region, period, mean_log2, n_samples)])
count_step("Gene rows (genes sharing transcript clusters combined)", nrow(rows))
count_step("Row x region x period cells", nrow(dev))
count_step("Cells with fewer than 3 samples (shown in grey)", sum(dev$n_samples < 3))

# developmental regulation (Kang et al. 2011 thresholds)
# expressed: mean log2 >= 6 in at least one region and period
# regulated: fetal (periods 3-7) vs postnatal (periods 8-12) brains. Each brain is the mean of its
#            samples across regions, so each brain counts once; Welch t-test on brain means;
#            BH FDR < 0.01 and |difference| >= 1 log2 (2-fold)
tst <- keep[period %in% c(3:12)]
reg <- rbindlist(lapply(seq_len(nrow(rows)), function(k) {
  ids <- intersect(ann[entrez == rows$entrez[k]]$ID, rownames(ex))
  d <- copy(tst)[, y := colMeans(ex[ids, sample, drop = FALSE], na.rm = TRUE)]
  b <- d[, .(y = mean(y), period = period[1]), by = donor]
  f <- b[period <= 7]$y; q <- b[period >= 8]$y
  data.table(row = rows$row[k], arm = rows$arm[k], n_fetal_brains = length(f), n_postnatal_brains = length(q),
             fetal_mean = mean(f), postnatal_mean = mean(q),
             log2_diff_fetal_vs_postnatal = mean(f) - mean(q),
             p = t.test(f, q)$p.value)
}))
mx <- dev[, .(max_mean_log2 = max(mean_log2)), by = row]
reg <- merge(reg, mx, by = "row")
reg[, fdr := p.adjust(p, "BH"), by = arm]   # mQTL and eQTL genes corrected separately
reg[, expressed := max_mean_log2 >= 6]
reg[, direction := fifelse(log2_diff_fetal_vs_postnatal > 0, "higher fetal", "higher postnatal")]
reg[, regulated := expressed & fdr < 0.01 & abs(log2_diff_fetal_vs_postnatal) >= 1]
fwrite(reg[order(-regulated, arm, row)], file.path(outdir, "dev_regulation.tsv"), sep = "\t")
for (a in c("mqtl", "eqtl")) {
  r <- reg[arm == a]
  count_step(sprintf("%s rows", a), nrow(r))
  count_step(sprintf("%s rows expressed (mean log2 >= 6 in a region and period)", a), sum(r$expressed))
  count_step(sprintf("%s rows developmentally regulated (expressed, FDR < 0.01, >= 2-fold)", a), sum(r$regulated))
  count_step(sprintf("%s   higher fetal", a), sum(r$regulated & r$direction == "higher fetal"))
  count_step(sprintf("%s   higher postnatal", a), sum(r$regulated & r$direction == "higher postnatal"))
}

# heatmaps: all periods, log2 expression
dev[, value := fifelse(n_samples >= 3, mean_log2, NA_real_)]
dev[, arm := factor(fifelse(arm == "mqtl", "mQTL", fifelse(arm == "eqtl", "eQTL", arm)), levels = c("mQTL", "eQTL"))]
dev[, region := factor(region, levels = c("NCX", "HIP", "AMY", "STR", "MD", "CBC"))]
dev[, row := factor(row, levels = rev(sort(unique(row))))]
short_label <- function(x) vapply(strsplit(as.character(x), "/"), function(g) {   # PCDHA1/.../PCDHA4 -> PCDHA1-4
  fam <- unique(sub("[0-9]+$", "", g)); num <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", g, perl = TRUE)))
  if (length(g) > 2 && length(fam) == 1 && !anyNA(num) && all(diff(sort(num)) == 1))
    sprintf("%s%d-%d", fam, min(num), max(num)) else paste(g, collapse = "/")
}, "")
plot_heatmap <- function(d, file) {
  p <- ggplot(d, aes(factor(period), row, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(. ~ region) +
    geom_vline(xintercept = 7.5, linewidth = 0.4, colour = "black", linetype = "dashed") +   # birth
    scale_x_discrete(breaks = as.character(seq(2, 14, 3))) +
    scale_y_discrete(labels = short_label) +
    scale_fill_viridis_c(name = expression(log[2] ~ "expression"), na.value = "grey80",
                         limits = range(dev$value, na.rm = TRUE)) +
    labs(x = "Developmental period (1-7 prenatal, 8-15 postnatal)", y = NULL) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 7, colour = "grey10"),
          strip.text = element_text(size = 8, face = "bold"), strip.text.y = element_text(angle = 0),
          legend.position = "bottom", legend.title = element_text(size = 7, vjust = 0.8),
          legend.text = element_text(size = 7), legend.key.height = unit(0.1, "in"),
          legend.key.width = unit(0.4, "in"), panel.spacing.x = unit(0.1, "in"))
  h <- 0.14 * uniqueN(d$row) + 1.0
  ggsave(paste0(file, ".pdf"), p, width = 6.7, height = h, device = pdf_device())
  ggsave(paste0(file, ".png"), p, width = 6.7, height = h, dpi = 300)
  saveRDS(p, paste0(file, ".rds"))
}
for (a in c("mQTL", "eQTL")) {
  d <- dev[arm == a]; d[, row := droplevels(row)]
  f <- file.path(outdir, paste0("dev_expression_heatmap_", tolower(a)))
  if (nrow(d)) plot_heatmap(d, f) else unlink(paste0(f, c(".pdf", ".png", ".rds")))
}
# trajectories of the developmentally regulated rows, one figure per dataset
for (a in c("mqtl", "eqtl")) {
  keep_rows <- reg[arm == a & regulated == TRUE][order(direction, row)]$row
  main <- dev[as.character(row) %in% keep_rows & !is.na(value)]
  if (!nrow(main)) {
    unlink(file.path(outdir, paste0("dev_trajectories_regulated_", a, c(".pdf", ".png", ".rds"))))
    next
  }
  main[, row := factor(as.character(row), levels = keep_rows)]
  p <- ggplot(main, aes(period, value, colour = region)) +
    geom_vline(xintercept = 7.5, linewidth = 0.3, linetype = "dashed") +   # birth
    geom_line(linewidth = 0.4) + geom_point(size = 0.7) +
    facet_wrap(~ row, ncol = 2) +
    scale_x_continuous(breaks = 1:15) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "Developmental period (1-7 prenatal, 8-15 postnatal)",
         y = expression(log[2] ~ "expression")) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 7, colour = "grey10"),
          strip.text = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
          legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  h <- 1.6 * ceiling(length(keep_rows) / 2) + 0.8
  ggsave(file.path(outdir, paste0("dev_trajectories_regulated_", a, ".pdf")), p, width = 6.7, height = h, device = pdf_device())
  ggsave(file.path(outdir, paste0("dev_trajectories_regulated_", a, ".png")), p, width = 6.7, height = h, dpi = 300)
  saveRDS(p, file.path(outdir, paste0("dev_trajectories_regulated_", a, ".rds")))
}
finish_log()
