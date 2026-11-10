# 16_flow_diagram.R
# Supplementary Figure S1: flow of the analysis from the input datasets to the prioritised
# features, with the methylation (mQTL) and expression (eQTL) arms side by side. All counts are
# read from the result files; the cytogenetic band is given when the prioritised features of an
# arm lie at one locus (org.Hs.eg.db MAP of the prioritised genes). Drawn at the final printed
# width (17 cm); text is 7-9 pt at that size.
#
# Output: results/figures/Supplementary_Figure_S1.pdf/.png
# Usage:  Rscript scripts/16_flow_diagram.R
# Requires ggplot2 and org.Hs.eg.db.

source("config/config.R")
source("scripts/functions.R")
suppressMessages({ library(ggplot2); library(org.Hs.eg.db) })
start_log("16_flow_diagram")
outdir <- file.path(PATHS$results, "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# counts
counts <- function(q) {
  a <- fread(file.path(PATHS$results, paste0("coloc_", q), "coloc_all_tests.tsv"))
  h <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits.tsv"))
  s <- fread(file.path(PATHS$results, paste0("smr_", q), "smr_results.tsv"))[analysis == "primary"]
  g <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits_gene_annotation.tsv"))
  pr <- s[smr_bonferroni_sig == TRUE & heidi_pass %in% TRUE, probeID]
  pr_loci <- unique(h[feature %in% pr, locus_id])
  band <- NA_character_
  if (length(pr_loci) == 1) {
    sym <- unlist(strsplit(na.omit(g[feature %in% pr, gene_symbol]), ";"))
    sym <- sym[nzchar(sym) & !grepl("^ENSG", sym)]
    m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, unique(sym), "MAP", "SYMBOL", multiVals = "first"))
    m <- m[!is.na(m)]
    if (length(m)) band <- names(sort(table(m), decreasing = TRUE))[1]
  }
  x <- list(loci_tested = uniqueN(a[status == "tested", locus_id]), tests = nrow(a[status == "tested"]),
            coloc = uniqueN(h$feature), coloc_loci = uniqueN(h$locus_id),
            testable = sum(!is.na(s$p_SMR)), sig = sum(s$smr_bonferroni_sig %in% TRUE),
            not_in_besd = sum(s$note %in% "probe not in BESD"),
            no_instrument = sum(s$note %in% "no instrument at this p threshold"),
            pass = length(pr), pr_loci = length(pr_loci), band = band)
  for (k in setdiff(names(x), "band")) count_step(sprintf("%s %s", q, k), x[[k]])
  logmsg(sprintf("  %s cytogenetic band of prioritised locus: %s", q, band))
  x
}
cm <- counts("mqtl")
ce <- counts("eqtl")

# drawing helpers (inches; y increases downwards)
pal <- list(gwas = c("#DCEAF7", "#5B8FC7"), ref = c("#ECECEC", "#8C8C8C"),
            m = c("#E6F3E4", "#5FA35B"), mF = c("#C8E6C3", "#3F8A3B"),
            e = c("#EEE6F6", "#8E6BBF"), eF = c("#DCCBEF", "#6E4AA3"),
            char = c("#FDEBD3", "#D9984A"))
pt <- function(x) x / .pt   # font size in points to ggplot size
lh <- 0.135; lh_side <- 0.125; pad <- 0.065
rects <- list(); texts <- list(); dashes <- list(); arrows <- list(); lines <- list()

box_height <- function(body, side) {
  2 * pad + lh * (1 + length(body)) + if (length(side)) 0.07 + lh_side * length(side) else 0
}
add_box <- function(x, y, w, h, colours, title, body = character(), side = character(),
                    side_parse = FALSE, lwd = 0.35) {
  rects[[length(rects) + 1]] <<- data.table(xmin = x - w / 2, xmax = x + w / 2, ymin = y, ymax = y + h,
                                            fill = colours[1], col = colours[2], lwd = lwd)
  n_main <- 1 + length(body)
  main_h <- lh * n_main + if (length(side)) 0.07 + lh_side * length(side) else 0
  y0 <- y + (h - main_h) / 2 + lh / 2
  texts[[length(texts) + 1]] <<- data.table(x = x, y = y0 + lh * (seq_len(n_main) - 1), label = c(title, body),
                                            face = c("bold", rep("plain", length(body))), size = c(pt(8), rep(pt(7.5), length(body))),
                                            parse = FALSE)
  if (length(side)) {
    ys <- y0 + lh * (n_main - 1) + lh / 2 + 0.03
    dashes[[length(dashes) + 1]] <<- data.table(x = x - w / 2 + 0.08, xend = x + w / 2 - 0.08, y = ys, col = colours[2])
    texts[[length(texts) + 1]] <<- data.table(x = x, y = ys + 0.04 + lh_side * (seq_along(side) - 0.5), label = side,
                                              face = "plain", size = pt(7), parse = side_parse)
  }
}
add_arrow <- function(x, y0, y1) arrows[[length(arrows) + 1]] <<- data.table(x = x, xend = x, y = y0, yend = y1)
add_line  <- function(x, y0, xend, y1) lines[[length(lines) + 1]] <<- data.table(x = x, xend = xend, y = y0, yend = y1)

W <- 6.7; bw <- 3.02; hgap <- 0.125; gap_y <- 0.15
cx <- W / 2; xm <- cx - (bw + hgap) / 2; xe <- cx + (bw + hgap) / 2

# input datasets
dw <- 1.53; dx <- cx + (-1.5:1.5) * (dw + hgap); y <- 0.04; dh <- box_height("x", character())
add_box(dx[1], y, dw, dh, pal$m, "Fetal brain mQTL", "Hannon et al. 2016")
add_box(dx[2], y, dw, dh, pal$gwas, "Schizophrenia GWAS", "PGC3, Trubetskoy et al. 2022")
add_box(dx[3], y, dw, dh, pal$e, "Fetal brain eQTL", "O'Brien et al. 2018")
add_box(dx[4], y, dw, dh, pal$ref, "1000 Genomes EUR", "Auton et al. 2015")
for (x in dx) add_arrow(x, y + dh, y + dh + 0.14)
y <- y + dh + 0.14

# quality control
qc_body <- c("GRCh37; alleles aligned to 1000 Genomes EUR, also the LD reference for fine-mapping and SMR",
             sprintf("Loci: lead SNP \u00b1 %d kb, overlapping windows merged, extended MHC excluded",
                     SETTINGS$locus_half_window / 1000))
qh <- box_height(qc_body, character())
add_box(cx, y, 4 * dw + 3 * hgap + 0.05, qh, pal$gwas, "Quality control and harmonisation", qc_body)
y <- y + qh
add_line(cx, y, cx, y + 0.08); add_line(xm, y + 0.08, xe, y + 0.08)
add_arrow(xm, y + 0.08, y + 0.2); add_arrow(xe, y + 0.08, y + 0.2)
y <- y + 0.2
arm_y <- y + 0.1
y <- y + 0.22

# one row of the two arms; boxes in a row share the height of the taller one
arm_row <- function(y, m, e, last = FALSE) {
  h <- max(box_height(m$body, m$side), box_height(e$body, e$side))
  add_box(xm, y, bw, h, m$pal, m$title, m$body, m$side, isTRUE(m$parse), if (is.null(m$lwd)) 0.35 else m$lwd)
  add_box(xe, y, bw, h, e$pal, e$title, e$body, e$side, isTRUE(e$parse), if (is.null(e$lwd)) 0.35 else e$lwd)
  if (!last) { add_arrow(xm, y + h, y + h + gap_y); add_arrow(xe, y + h, y + h + gap_y) }
  y + h + if (last) 0 else gap_y
}
coloc_txt <- "coloc.susie (coloc.abf if no credible set)"
y <- arm_row(y, list(pal = pal$m, title = "Colocalisation per CpG", body = c(coloc_txt, sprintf("%d loci, %d CpG tests", cm$loci_tested, cm$tests))),
                list(pal = pal$e, title = "Colocalisation per gene", body = c(coloc_txt, sprintf("%d loci, %d gene tests", ce$loci_tested, ce$tests))))
pp4_title <- sprintf("PP4 \u2265 %s", SETTINGS$pp4_threshold)
y <- arm_row(y, list(pal = pal$m, title = pp4_title, body = sprintf("%d CpGs colocalised at %d loci", cm$coloc, cm$coloc_loci)),
                list(pal = pal$e, title = pp4_title, body = sprintf("%d genes colocalised at %d loci", ce$coloc, ce$coloc_loci)))

not_testable <- function(x, unit, qtl) {
  thr <- strsplit(formatC(SETTINGS$smr_peqtl_main, format = "e", digits = 0), "e")[[1]]
  out <- character(); parse <- FALSE
  if (x$not_in_besd > 0)
    out <- c(out, sprintf("Not testable: %d %s absent from the SMR-format %s data", x$not_in_besd, unit, qtl))
  if (x$no_instrument > 0) {
    out <- c(out, sprintf("Not~testable*':'~%d~%s~with~no~%s~SNP~at~p < %s %%*%% 10^%d",
                          x$no_instrument, unit, qtl, as.integer(thr[1]), as.integer(thr[2])))
    parse <- TRUE
  }
  if (parse && x$not_in_besd > 0) stop("Both reasons in one arm: adjust the side text")
  list(text = out, parse = parse)
}
ntm <- not_testable(cm, "CpGs", "mQTL"); nte <- not_testable(ce, "genes", "eQTL")
y <- arm_row(y, list(pal = pal$m, title = "SMR, Bonferroni-significant", body = sprintf("%d of %d testable CpGs", cm$sig, cm$testable),
                     side = ntm$text, parse = ntm$parse),
                list(pal = pal$e, title = "SMR, Bonferroni-significant", body = sprintf("%d of %d testable genes", ce$sig, ce$testable),
                     side = nte$text, parse = nte$parse))
y <- arm_row(y, list(pal = pal$m, title = "HEIDI p > 0.05", body = sprintf("%d pass", cm$pass)),
                list(pal = pal$e, title = "HEIDI p > 0.05", body = sprintf("%d pass", ce$pass)))
prioritised <- function(x, unit, table) {
  if (x$pass == 0) return(list(title = sprintf("Prioritised: 0 %s", unit), body = character()))
  at <- if (x$pr_loci == 1 && !is.na(x$band)) x$band else sprintf("%d loci", x$pr_loci)
  list(title = sprintf("Prioritised: %d %s at %s", x$pass, unit, at), body = sprintf("(%s)", table))
}
pm <- prioritised(cm, "CpGs", "Table 2"); pe <- prioritised(ce, "genes", "Table 3")
y <- arm_row(y, list(pal = pal$mF, title = pm$title, body = pm$body, lwd = 0.55),
                list(pal = pal$eF, title = pe$title, body = pe$body, lwd = 0.55),
             last = TRUE)

# characterisation
add_line(xm, y, xm, y + 0.08); add_line(xe, y, xe, y + 0.08); add_line(xm, y + 0.08, xe, y + 0.08)
add_arrow(cx, y + 0.08, y + 0.2)
y <- y + 0.2
ch_body <- c("Developmental expression (Kang et al. 2011)", "Cell-type expression (Polioudakis et al. 2019)")
chh <- box_height(ch_body, character())
add_box(cx, y, 4.8, chh, pal$char, "Characterisation of prioritised genes", ch_body)
H <- y + chh + 0.08

# plot
R <- rbindlist(rects); Tx <- rbindlist(texts); D <- rbindlist(dashes); A <- rbindlist(arrows); L <- rbindlist(lines)
arms <- data.table(x = c(xm, xe), y = arm_y, label = c("DNA methylation (mQTL)", "Gene expression (eQTL)"),
                   col = c("#2F6E2C", "#5A3A8C"))
p <- ggplot() +
  geom_rect(data = R, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin, fill = fill, colour = col, linewidth = lwd)) +
  geom_segment(data = L, aes(x = x, xend = xend, y = -y, yend = -yend), colour = "grey20", linewidth = 0.35) +
  geom_segment(data = A, aes(x = x, xend = xend, y = -y, yend = -yend), colour = "grey20", linewidth = 0.35,
               arrow = arrow(length = unit(0.055, "in"), type = "closed")) +
  geom_text(data = Tx[parse == FALSE], aes(x, -y, label = label, fontface = face, size = size), colour = "grey10") +
  geom_text(data = Tx[parse == TRUE], aes(x, -y, label = label, size = size), colour = "grey10", parse = TRUE) +
  geom_text(data = arms, aes(x, -y, label = label, colour = col), fontface = "bold", size = pt(9)) +
  scale_fill_identity() + scale_colour_identity() + scale_linewidth_identity() + scale_size_identity() +
  coord_fixed(xlim = c(0, W), ylim = c(-H, 0), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
if (nrow(D)) p <- p + geom_segment(data = D, aes(x = x, xend = xend, y = -y, yend = -y, colour = col),
                                   linetype = "dashed", linewidth = 0.25)

ggsave(file.path(outdir, "Supplementary_Figure_S1.pdf"), p, width = W, height = H, device = pdf_device())
ggsave(file.path(outdir, "Supplementary_Figure_S1.png"), p, width = W, height = H, dpi = 300, bg = "white")
logmsg("  Figure written to ", file.path(outdir, "Supplementary_Figure_S1.pdf"))
finish_log()
