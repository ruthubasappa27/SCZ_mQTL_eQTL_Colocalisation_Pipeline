# 03_qc_mqtl.R
# Quality control of the Hannon et al. (2016) fetal brain mQTL data and alignment to the
# 1000 Genomes EUR reference panel. Every SNP-CpG pair is kept.
#
# Output: data/processed/mqtl_harmonised.tsv.gz  (FEATURE = CpG; BETA = effect of panel allele A1)
# Usage:  Rscript scripts/03_qc_mqtl.R
#
# Notes:
#   - SNP_Allele is the allele the beta refers to; the other allele is the second panel allele
#     at the same position.
#   - Standard errors are not provided and are derived from beta and the two-sided p-value.

source("config/config.R")
source("scripts/functions.R")
start_log("03_qc_mqtl")

m <- fread(PATHS$mqtl_raw)
need_cols(m, c("ProbeID", "SNP_Chr", "SNP_BP", "SNP_Allele", "beta", "p.value"), "mQTL file")
count_step("mQTL: raw SNP-CpG pairs", nrow(m))
count_step("mQTL: raw unique CpGs", uniqueN(m$ProbeID))

opt <- function(col) if (col %in% names(m)) m[[col]] else NA
m <- data.table(FEATURE = m$ProbeID,
                FEATURE_CHR = opt("DNAm_CHR"), FEATURE_POS = opt("DNAm_BP"),
                FEATURE_LABEL = opt("GeneAnnotation"),
                CHR = as.integer(m$SNP_Chr), POS = as.integer(m$SNP_BP),
                EA = toupper(m$SNP_Allele), BETA_raw = m$beta, P = m$p.value)

m <- m[is.finite(BETA_raw) & is.finite(P)];          count_step("mQTL: non-missing beta and p", nrow(m))
count_step("mQTL: p-value exactly 0 (SE cannot be derived; removed)", sum(m$P == 0))
m <- m[P > 0 & P <= 1 & BETA_raw != 0]

# SE = |beta| / Z, with Z from the two-sided p-value
m[, SE := abs(BETA_raw) / qnorm(P / 2, lower.tail = FALSE)]
m <- m[is.finite(SE) & SE > 0];                      count_step("mQTL: valid derived SE", nrow(m))

# match to the reference panel by position and allele
bim <- load_bim()
mm <- merge(m, bim[, .(CHR, POS, SNP, A1, A2)], by = c("CHR", "POS"), allow.cartesian = TRUE)
mm <- mm[EA == A1 | EA == A2]
mm[, n_map := .N, by = .(FEATURE, CHR, POS)]
count_step("mQTL: pairs matching a panel SNP by position + allele", uniqueN(mm[, .(FEATURE, CHR, POS)]))
count_step("mQTL: ambiguous (>1 panel SNP at position; removed)", uniqueN(mm[n_map > 1, .(FEATURE, CHR, POS)]))
mm <- mm[n_map == 1]

mm[, BETA := fifelse(EA == A1, BETA_raw, -BETA_raw)]
count_step("mQTL: effect allele = panel A2 (beta sign flipped)", sum(mm$EA == mm$A2))

# duplicate records of one SNP-CpG pair: keep the most significant
setorder(mm, FEATURE, SNP, P)
n0 <- nrow(mm)
mm <- mm[!duplicated(mm[, .(FEATURE, SNP)])]
count_step("mQTL: duplicate SNP-CpG records removed", n0 - nrow(mm))
count_step("mQTL: final SNP-CpG pairs", nrow(mm))
count_step("mQTL: final unique CpGs", uniqueN(mm$FEATURE))
count_step("mQTL: final unique SNPs", uniqueN(mm$SNP))

fwrite(mm[, .(FEATURE, FEATURE_LABEL, FEATURE_CHR, FEATURE_POS,
              SNP, CHR, POS, A1, A2, EA_original = EA, BETA, SE, P)],
       file.path(PATHS$processed, "mqtl_harmonised.tsv.gz"), sep = "\t")
finish_log()
