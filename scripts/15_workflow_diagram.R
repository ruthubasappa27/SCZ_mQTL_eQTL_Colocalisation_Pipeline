# 15_workflow_diagram.R
# Figure 1: analysis workflow (methods only, no counts). Colour marks the datasets only
# (GWAS, mQTL, eQTL, reference panel); shared steps are neutral. Drawn at the final printed width
# (17 cm); text is 7.5-9 pt at that size. The version with counts
# (Supplementary Figure S1) is drawn by 16_flow_diagram.R.
#
# Output: results/figures/Figure1.pdf/.png
# Usage:  Rscript scripts/15_workflow_diagram.R
# Requires ggplot2.

source("config/config.R")
source("scripts/functions.R")
suppressMessages(library(ggplot2))
start_log("15_workflow_diagram")
outdir <- file.path(PATHS$results, "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# fill and border per dataset; shared steps use one colour
pal <- list(gwas = c("#DCEAF7", "#5B8FC7"), ref = c("#ECECEC", "#8C8C8C"),
            m = c("#E6F3E4", "#5FA35B"), e = c("#EEE6F6", "#8E6BBF"),
            step = c("#F6F4EF", "#6B6B6B"), final = c("#F6F4EF", "#333333"))
arm_text <- c("#2F6E2C", "#5A3A8C")

pt <- function(x) x / .pt
lh <- 0.135; pad <- 0.065
rects <- list(); texts <- list(); arrows <- list(); lines <- list()

box_height <- function(body) 2 * pad + lh * (1 + length(body))
add_box <- function(x, y, w, colours, title, body = character(), lwd = 0.35) {
  h <- box_height(body)
  rects[[length(rects) + 1]] <<- data.table(xmin = x - w / 2, xmax = x + w / 2, ymin = y, ymax = y + h,
                                            fill = colours[1], col = colours[2], lwd = lwd)
  n <- 1 + length(body)
  texts[[length(texts) + 1]] <<- data.table(x = x, y = y + pad + lh * (seq_len(n) - 0.5), label = c(title, body),
                                            face = c("bold", rep("plain", length(body))),
                                            size = c(pt(8), rep(pt(7.5), length(body))))
  y + h
}
add_arrow <- function(x, y0, y1) arrows[[length(arrows) + 1]] <<- data.table(x = x, xend = x, y = y0, yend = y1)
add_line  <- function(x, y0, xend, y1) lines[[length(lines) + 1]] <<- data.table(x = x, xend = xend, y = y0, yend = y1)

W <- 6.7; bw <- 3.02; hgap <- 0.125; cw <- 4.8; ah <- 0.14
cx <- W / 2; xm <- cx - (bw + hgap) / 2; xe <- cx + (bw + hgap) / 2

# input datasets
dw <- 1.53; dx <- cx + (-1.5:1.5) * (dw + hgap); y <- 0.04
add_box(dx[1], y, dw, pal$m, "Fetal brain mQTL", "Hannon et al. 2016")
add_box(dx[2], y, dw, pal$gwas, "Schizophrenia GWAS", "PGC3, Trubetskoy et al. 2022")
add_box(dx[3], y, dw, pal$e, "Fetal brain eQTL", "O'Brien et al. 2018")
y <- add_box(dx[4], y, dw, pal$ref, "1000 Genomes EUR", "Auton et al. 2015")
for (x in dx) add_arrow(x, y, y + ah)
y <- y + ah

# quality control
y <- add_box(cx, y, 4 * dw + 3 * hgap + 0.05, pal$step, "Quality control and harmonisation",
             c("GRCh37; alleles aligned to 1000 Genomes EUR, also the LD reference for fine-mapping and SMR",
               sprintf("Loci: lead SNP \u00b1 %d kb, overlapping windows merged, extended MHC excluded",
                       SETTINGS$locus_half_window / 1000)))
add_line(cx, y, cx, y + 0.08); add_line(xm, y + 0.08, xe, y + 0.08)
add_arrow(xm, y + 0.08, y + 0.2); add_arrow(xe, y + 0.08, y + 0.2)
arm_y <- y + 0.3
y <- y + 0.42

# colocalisation, one box per dataset
coloc_body <- c("coloc.susie (coloc.abf if no credible set)", sprintf("PP4 \u2265 %s", SETTINGS$pp4_threshold))
add_box(xm, y, bw, pal$m, "Colocalisation per CpG", coloc_body)
y <- add_box(xe, y, bw, pal$e, "Colocalisation per gene", coloc_body)
add_line(xm, y, xm, y + 0.08); add_line(xe, y, xe, y + 0.08); add_line(xm, y + 0.08, xe, y + 0.08)
add_arrow(cx, y + 0.08, y + 0.2)
y <- y + 0.2

# shared steps
y <- add_box(cx, y, cw, pal$step, "Summary-data-based Mendelian randomisation",
             "SMR, Bonferroni-corrected; HEIDI p > 0.05")
add_arrow(cx, y, y + ah); y <- y + ah
y <- add_box(cx, y, cw, pal$final, "Prioritised CpGs and genes",
             "colocalised, SMR-significant and passing HEIDI", lwd = 0.6)
add_arrow(cx, y, y + ah); y <- y + ah
y <- add_box(cx, y, cw, pal$step, "Characterisation of prioritised genes",
             c("Developmental expression (Kang et al. 2011)", "Cell-type expression (Polioudakis et al. 2019)"))
H <- y + 0.08

# plot
R <- rbindlist(rects); Tx <- rbindlist(texts); A <- rbindlist(arrows); L <- rbindlist(lines)
arms <- data.table(x = c(xm, xe), y = arm_y, label = c("DNA methylation (mQTL)", "Gene expression (eQTL)"), col = arm_text)
p <- ggplot() +
  geom_rect(data = R, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin, fill = fill, colour = col, linewidth = lwd)) +
  geom_segment(data = L, aes(x = x, xend = xend, y = -y, yend = -yend), colour = "grey20", linewidth = 0.35) +
  geom_segment(data = A, aes(x = x, xend = xend, y = -y, yend = -yend), colour = "grey20", linewidth = 0.35,
               arrow = arrow(length = unit(0.055, "in"), type = "closed")) +
  geom_text(data = Tx, aes(x, -y, label = label, fontface = face, size = size), colour = "grey10") +
  geom_text(data = arms, aes(x, -y, label = label, colour = col), fontface = "bold", size = pt(9)) +
  scale_fill_identity() + scale_colour_identity() + scale_linewidth_identity() + scale_size_identity() +
  coord_fixed(xlim = c(0, W), ylim = c(-H, 0), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(outdir, "Figure1.pdf"), p, width = W, height = H, device = pdf_device())
ggsave(file.path(outdir, "Figure1.png"), p, width = W, height = H, dpi = 300, bg = "white")
logmsg("  Figure written to ", file.path(outdir, "Figure1.pdf"))
finish_log()
