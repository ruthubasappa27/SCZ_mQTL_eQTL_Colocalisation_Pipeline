# 17_main_figures.R
# Figures 2 and 3 of the manuscript, assembled from the plots saved by 12 and 13 (mQTL genes).
#   Figure 2: A, developmental expression heatmap; B, trajectories of the developmentally
#             regulated genes.
#   Figure 3: cell-type expression heatmap (single panel).
# Output: results/figures/Figure2.pdf/.png, Figure3.pdf/.png
# Usage:  Rscript scripts/17_main_figures.R

source("config/config.R")
source("scripts/functions.R")
suppressMessages({ library(ggplot2); library(grid); library(gtable) })
start_log("17_main_figures")
indir  <- file.path(PATHS$results, "characterisation")
outdir <- file.path(PATHS$results, "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_plot <- function(name) {
  f <- file.path(indir, paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing ", f, "; run 12 and 13 first.")
  readRDS(f)
}

panel_label <- function(x)
  textGrob(x, x = unit(0, "npc"), hjust = 0, gp = gpar(fontsize = 10, fontface = "bold"))

# heights in inches; widths as a fraction of the figure width (panel left-aligned);
# panel letters only when a figure has more than one panel
assemble <- function(plots, heights, widths, file, width = 6.7) {
  lab_h <- if (length(plots) > 1) 0.2 else 0
  g <- gtable(widths = unit(1, "null"), heights = unit(as.vector(rbind(lab_h, heights)), "in"))
  for (i in seq_along(plots)) {
    p <- ggplotGrob(plots[[i]])
    if (widths[i] < 1) p <- gTree(children = gList(p), vp = viewport(x = 0, width = widths[i], just = "left"))
    if (lab_h > 0)
      g <- gtable_add_grob(g, panel_label(LETTERS[i]), t = 2 * i - 1, l = 1, name = paste0("label", i))
    g <- gtable_add_grob(g, p, t = 2 * i, l = 1, name = paste0("panel", i))
  }
  h <- sum(heights) + lab_h * length(plots)
  ggsave(paste0(file, ".pdf"), g, width = width, height = h, device = pdf_device())
  ggsave(paste0(file, ".png"), g, width = width, height = h, dpi = 300, bg = "white")
  logmsg("  ", basename(file), ": ", width, " x ", round(h, 2), " in")
}

n_rows <- function(p, col) length(unique(p$data[[col]]))

# Figure 2
hm  <- read_plot("dev_expression_heatmap_mqtl")
trj <- read_plot("dev_trajectories_regulated_mqtl")
assemble(list(hm, trj),
         heights = c(0.14 * n_rows(hm, "row") + 1.0, 1.6 * ceiling(n_rows(trj, "row") / 2) + 0.8),
         widths  = c(1, 1),
         file    = file.path(outdir, "Figure2"))

# Figure 3
ct <- read_plot("celltype_heatmap_mqtl")
assemble(list(ct),
         heights = 0.15 * n_rows(ct, "Gene") + 0.9,
         widths  = 1,
         file    = file.path(outdir, "Figure3"))

finish_log()
