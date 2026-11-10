# 01_qc_gwas.R
# Quality control of the PGC3 schizophrenia GWAS summary statistics and alignment to the
# 1000 Genomes EUR reference panel.
#
# Output: data/processed/gwas_qc.tsv.gz          QC'd variants (used to define loci)
#         data/processed/gwas_harmonised.tsv.gz  aligned to the panel; BETA = effect of panel allele A1
# Usage:  Rscript scripts/01_qc_gwas.R

source("config/config.R")
source("scripts/functions.R")
start_log("01_qc_gwas")
dir.create(PATHS$processed, showWarnings = FALSE, recursive = TRUE)

g <- fread(PATHS$gwas_raw, skip = "CHROM")
setnames(g, sub("^#", "", names(g)))
neff_col <- intersect(c("NEFF", "NEFFDIV2"), names(g))[1]
need_cols(g, c("CHROM", "ID", "POS", "A1", "A2", "FCAS", "FCON", "IMPINFO",
               "BETA", "SE", "PVAL", "NCAS", "NCON"), "GWAS file")
if (is.na(neff_col)) stop("GWAS file has no NEFF/NEFFDIV2 column.")
count_step("GWAS: raw variants", nrow(g))

# N = total effective sample size (NEFF in the PGC3 header; NEFFDIV2 is half of it)
g <- g[, .(CHR = as.integer(CHROM), SNP = ID, POS = as.integer(POS),
           EA = toupper(A1), OA = toupper(A2),
           FCAS, FCON, NCAS, NCON, INFO = IMPINFO,
           BETA, SE, P = PVAL,
           N = get(neff_col) * ifelse(neff_col == "NEFFDIV2", 2, SETTINGS$gwas_neff_multiplier))]
logmsg("  Effective total N = ", neff_col, " x ", SETTINGS$gwas_neff_multiplier,
       " (median ", round(median(g$N)), ")")

g <- g[CHR %in% 1:22];                               count_step("GWAS: autosomal", nrow(g))
g <- g[INFO >= SETTINGS$gwas_info_min];              count_step("GWAS: INFO >= 0.8", nrow(g))

# effect-allele frequency in cases and controls combined
g[, EAF := (FCAS * NCAS + FCON * NCON) / (NCAS + NCON)]
g[, MAF := pmin(EAF, 1 - EAF)]
g <- g[MAF >= SETTINGS$maf_min];                     count_step("GWAS: MAF >= 0.01 (cases+controls)", nrow(g))
g <- g[is_snv(EA) & is_snv(OA)];                     count_step("GWAS: biallelic SNVs (A/C/G/T)", nrow(g))
if (SETTINGS$drop_palindromic) {
  g <- g[!is_palindromic(EA, OA)];                   count_step("GWAS: palindromic A/T, C/G removed", nrow(g))
}
g <- g[is.finite(SE) & SE > 0 & is.finite(BETA) & abs(BETA) < SETTINGS$max_abs_beta &
       is.finite(P) & P > 0 & P <= 1];               count_step("GWAS: valid BETA/SE/P", nrow(g))
dup <- g$SNP[duplicated(g$SNP)]
g <- g[!SNP %in% dup];                               count_step("GWAS: duplicated IDs removed (all copies)", nrow(g))

fwrite(g[, .(SNP, CHR, POS, EA, OA, EAF, MAF, BETA, SE, P, N)],
       file.path(PATHS$processed, "gwas_qc.tsv.gz"), sep = "\t")

# alignment to the reference panel
bim <- load_bim()
h <- merge(g, bim[, .(SNP, CHR_ref = CHR, POS_ref = POS, A1, A2)], by = "SNP")
count_step("GWAS: rsID present in 1000G EUR panel", nrow(h))
h <- h[CHR == CHR_ref & POS == POS_ref];             count_step("GWAS: same chromosome/position as panel", nrow(h))
h[, align := fifelse(EA == A1 & OA == A2, "same",
             fifelse(EA == A2 & OA == A1, "swapped", "mismatch"))]
count_step("GWAS: alleles do not match panel (removed)", sum(h$align == "mismatch"))
h <- h[align != "mismatch"]
count_step("GWAS: alleles swapped to panel A1 (beta sign flipped)", sum(h$align == "swapped"))
h[align == "swapped", `:=`(BETA = -BETA, EAF = 1 - EAF)]
count_step("GWAS: final harmonised SNPs", nrow(h))

fwrite(h[, .(SNP, CHR, POS, A1, A2, FREQ_A1 = EAF, MAF, BETA, SE, P, N)],
       file.path(PATHS$processed, "gwas_harmonised.tsv.gz"), sep = "\t")
finish_log()
