# 02_define_loci.R
# Defines GWAS loci: lead SNPs (p < 5e-8) are selected greedily by p-value, removing all SNPs
# within 500 kb of each lead; each lead defines a +/- 500 kb window; overlapping windows are
# merged. Windows overlapping the extended MHC (chr6:25-34 Mb, GRCh37) form one locus, which
# is excluded or kept according to SETTINGS$mhc_action.
#
# Output: results/loci/lead_snps.tsv, results/loci/loci_all.tsv, results/loci/loci_for_coloc.tsv
# Usage:  Rscript scripts/02_define_loci.R

source("config/config.R")
source("scripts/functions.R")
start_log("02_define_loci")
outdir <- file.path(PATHS$results, "loci")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
hw <- SETTINGS$locus_half_window

g <- fread(file.path(PATHS$processed, "gwas_qc.tsv.gz"))
sig <- g[P < SETTINGS$gwas_p_sig][order(P)]
count_step("Genome-wide significant SNPs", nrow(sig))

leads <- list()
while (nrow(sig) > 0) {
  top <- sig[1]
  leads[[length(leads) + 1]] <- top
  sig <- sig[!(CHR == top$CHR & abs(POS - top$POS) <= hw)]
}
leads <- rbindlist(leads)
count_step("Lead SNPs (greedy, 500 kb)", nrow(leads))
fwrite(leads, file.path(outdir, "lead_snps.tsv"), sep = "\t")

w <- leads[, .(CHR, start = pmax(1L, POS - hw), end = POS + hw, lead_snp = SNP, lead_p = P)]
w[, in_mhc := CHR == SETTINGS$mhc_chr & end >= SETTINGS$mhc_start & start <= SETTINGS$mhc_end]
count_step("Lead-SNP windows touching the extended MHC", sum(w$in_mhc))

# merge overlapping windows; MHC windows are merged into one block separately
merge_windows <- function(x) {
  x <- copy(x)
  setorder(x, CHR, start)
  x[, run_end := cummax(end), by = CHR]
  x[, new_block := c(TRUE, start[-1] > run_end[-.N]), by = CHR]
  x[, grp := cumsum(new_block)]
  x[, .(CHR = CHR[1], start = min(start), end = max(end),
        lead_snp = lead_snp[which.min(lead_p)], lead_p = min(lead_p),
        n_lead_snps = .N, merged_leads = paste(lead_snp, collapse = ",")), by = grp][, grp := NULL][]
}
loci <- merge_windows(w[in_mhc == FALSE])
loci[, in_mhc := FALSE]
if (any(w$in_mhc)) {
  mhc <- w[in_mhc == TRUE][, .(CHR = CHR[1], start = min(start), end = max(end),
                               lead_snp = lead_snp[which.min(lead_p)], lead_p = min(lead_p),
                               n_lead_snps = .N, merged_leads = paste(lead_snp, collapse = ","))]
  mhc[, in_mhc := TRUE]
  loci <- rbind(loci, mhc)
}
setorder(loci, CHR, start)
loci[, locus_id := sprintf("chr%d_%d_%d", CHR, start, end)]
loci[, status := fifelse(in_mhc & SETTINGS$mhc_action == "exclude", "excluded_MHC", "analysed")]
count_step("Loci after merging overlapping windows (incl. MHC block)", nrow(loci))
count_step("Loci formed by merging >1 lead-SNP window", sum(loci$n_lead_snps > 1))
count_step("Loci taken forward to colocalisation", sum(loci$status == "analysed"))

setcolorder(loci, c("locus_id", "CHR", "start", "end", "lead_snp", "lead_p",
                    "n_lead_snps", "in_mhc", "status", "merged_leads"))
fwrite(loci, file.path(outdir, "loci_all.tsv"), sep = "\t")
fwrite(loci[status == "analysed"], file.path(outdir, "loci_for_coloc.tsv"), sep = "\t")
finish_log()
