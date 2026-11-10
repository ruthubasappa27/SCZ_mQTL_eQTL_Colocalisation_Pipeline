# 13_cell_type_expression.R
# Cell-type expression of the prioritised genes in the mid-gestation human neocortex atlas of
# Polioudakis et al. (2019), from CoDEx summary statistics (mean normalised expression per cell
# type), reported separately for the mQTL and eQTL genes.
#   Input:  data/codex/codex_mqtl.csv and data/codex/codex_eqtl.csv, exported from the CoDEx
#           heatmap summary (columns Gene, Mean, SD, Percent, Cluster_name; Cluster_name NA = all
#           cells). Genes are matched to the prioritised genes by current symbol
#           (results/characterisation/dev_expression_gene_availability.tsv); genes CoDEx did not
#           recognise are reported as not available.
#   Per gene: cell type with the highest mean, that mean, and the range (highest minus lowest
#           mean across the 16 cell types). Range plots of both datasets share one x axis.
# Figures are drawn at printed size (heatmaps 17 cm, range plots 8.5 cm wide; text >= 7 pt).
# Plot objects are saved (.rds) for 17_main_figures.R.
#
# Output: results/characterisation/celltype_summary.tsv, celltype_values.tsv,
#         celltype_heatmap_<qtl>.pdf/.png/.rds (all genes), celltype_range_<qtl>.pdf/.png/.rds
# Usage:  Rscript scripts/13_cell_type_expression.R

source("config/config.R")
source("scripts/functions.R")
suppressMessages(library(ggplot2))
start_log("13_cell_type_expression")
outdir <- file.path(PATHS$results, "characterisation")
indir  <- PATHS$codex_dir

cell_types <- c("vRG" = "vRG", "oRG" = "oRG",
                "Cycling progenitors (S phase)" = "PgS", "Cycling progenitors (G2/M phase)" = "PgG2M",
                "IP" = "IP", "Migrating excitatory" = "ExN", "Maturing excitatory" = "ExM",
                "Maturing excitatory upper enriched" = "ExM-U", "Excitatory deep layer 1" = "ExDp1",
                "Excitatory Deep layer 2" = "ExDp2", "Interneuron MGE" = "InMGE",
                "Interneuron CGE" = "InCGE", "OPCs" = "OPC", "Endothelial" = "End",
                "Pericyte" = "Per", "Microglia" = "Mic")

av <- fread(file.path(outdir, "dev_expression_gene_availability.tsv"))
av[is.na(current_symbol) | current_symbol == "", current_symbol := gene]

arm_colour <- c(mqtl = "#3F8A3B", eqtl = "#6E4AA3")

arms <- intersect(c("mqtl", "eqtl"), unique(av$qtl))
for (q in setdiff(c("mqtl", "eqtl"), arms)) {
  logmsg("  ", q, ": no prioritised genes; no cell-type analysis")
  unlink(file.path(outdir, paste0(rep(c("celltype_heatmap_", "celltype_range_"), 3), q,
                                  rep(c(".pdf", ".png", ".rds"), each = 2))))
}

# common x axis and colour scale for both datasets
all_means <- unlist(lapply(arms, function(q) {
  f <- file.path(indir, paste0("codex_", q, ".csv")); if (file.exists(f)) fread(f)$Mean }))
xlim <- range(all_means, na.rm = TRUE)

summ <- list(); vals <- list()
for (q in arms) {
  f <- file.path(indir, paste0("codex_", q, ".csv"))
  if (!file.exists(f)) stop("Missing CoDEx export: ", f)
  d <- fread(f)
  need_cols(d, c("Gene", "Mean", "SD", "Percent", "Cluster_name"), basename(f))
  d <- d[!is.na(Cluster_name) & Cluster_name != "NA"]
  unknown <- setdiff(unique(d$Cluster_name), names(cell_types))
  if (length(unknown)) stop("Unrecognised CoDEx cell types in ", basename(f), ": ", paste(unknown, collapse = ", "))
  d[, cell_type := factor(cell_types[Cluster_name], levels = cell_types)]

  genes_q <- unique(av[qtl == q, .(gene, current_symbol)])
  extra <- setdiff(unique(d$Gene), genes_q$current_symbol)
  if (length(extra)) stop("Genes in ", basename(f), " that are not prioritised ", q, " genes: ", paste(extra, collapse = ", "))
  bad <- d[, .N, by = Gene][N != length(cell_types)]
  if (nrow(bad)) stop("Genes without all 16 cell types in ", basename(f), ": ", paste(bad$Gene, collapse = ", "))

  s <- d[, .(peak_cell_type = as.character(cell_type[which.max(Mean)]), peak_mean = max(Mean),
             min_mean = min(Mean), range = max(Mean) - min(Mean)), by = .(current_symbol = Gene)]
  s <- merge(genes_q, s, by = "current_symbol", all.x = TRUE)
  s[, `:=`(qtl = q, in_codex = !is.na(peak_mean))]
  count_step(sprintf("%s prioritised genes", q), nrow(s))
  count_step(sprintf("%s genes in CoDEx", q), sum(s$in_codex))
  count_step(sprintf("%s genes absent from the CoDEx export", q), sum(!s$in_codex))
  if (any(!s$in_codex)) logmsg("    absent: ", paste(s[in_codex == FALSE, current_symbol], collapse = ", "))
  summ[[q]] <- s
  vals[[q]] <- cbind(qtl = q, d[, .(gene = Gene, cell_type, mean = Mean, sd = SD, percent = Percent)])

  # heatmap of all genes (mean normalised expression)
  d[, Gene := factor(Gene, levels = rev(sort(unique(Gene))))]
  p <- ggplot(d, aes(cell_type, Gene, fill = Mean)) + geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(name = "Mean normalised\nexpression", limits = xlim) +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 7, colour = "grey10"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          legend.key.width = unit(0.12, "in"))
  h <- max(2, 0.15 * uniqueN(d$Gene) + 0.9)
  ggsave(file.path(outdir, paste0("celltype_heatmap_", q, ".pdf")), p, width = 6.7, height = h, device = pdf_device())
  ggsave(file.path(outdir, paste0("celltype_heatmap_", q, ".png")), p, width = 6.7, height = h, dpi = 300)
  saveRDS(p, file.path(outdir, paste0("celltype_heatmap_", q, ".rds")))

  # range plot (compact): lowest to highest cell-type mean per gene, peak cell type labelled
  r <- s[in_codex == TRUE][order(peak_mean)]
  r[, current_symbol := factor(current_symbol, levels = current_symbol)]
  p <- ggplot(r, aes(y = current_symbol)) +
    geom_segment(aes(x = min_mean, xend = peak_mean, yend = current_symbol), colour = "grey60", linewidth = 0.4) +
    geom_point(aes(x = min_mean), colour = "grey40", size = 0.9) +
    geom_point(aes(x = peak_mean), colour = arm_colour[[q]], size = 1.3) +
    geom_text(aes(x = peak_mean, label = peak_cell_type), hjust = -0.25, size = 7 / .pt) +
    scale_x_continuous(limits = c(xlim[1], xlim[2] * 1.2)) +
    labs(x = "Mean normalised expression", y = NULL) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 7, colour = "grey10"))
  h <- 0.14 * nrow(r) + 0.6
  ggsave(file.path(outdir, paste0("celltype_range_", q, ".pdf")), p, width = 3.35, height = h, device = pdf_device())
  ggsave(file.path(outdir, paste0("celltype_range_", q, ".png")), p, width = 3.35, height = h, dpi = 300)
  saveRDS(p, file.path(outdir, paste0("celltype_range_", q, ".rds")))
}
summ <- rbindlist(summ, fill = TRUE)
fwrite(summ[order(qtl, -peak_mean)], file.path(outdir, "celltype_summary.tsv"), sep = "\t")
fwrite(rbindlist(vals), file.path(outdir, "celltype_values.tsv"), sep = "\t")
print(summ[order(qtl, -peak_mean), .(qtl, gene = current_symbol, peak_cell_type,
                                     peak_mean = round(peak_mean, 3), range = round(range, 3))])
finish_log()
