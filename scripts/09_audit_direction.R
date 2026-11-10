# 09_audit_direction.R
# Re-derives the direction of effect of the risk allele at the top SNP of every colocalised
# feature from the original (unharmonised) GWAS, mQTL and eQTL files, and compares it with
# the result of 06_results.R.
#   GWAS: BETA = log(OR) of A1
#   mQTL: beta refers to SNP_Allele
#   eQTL: slope per ALT allele; REF/ALT from snp_positions.bed
#
# Output: results/direction_audit.tsv
# Usage:  Rscript scripts/09_audit_direction.R (after 04_qc_eqtl.R, which writes the eQTL subset files)

source("config/config.R")
source("scripts/functions.R")
start_log("09_audit_direction")

hm <- fread(file.path(PATHS$results, "coloc_mqtl", "hits_directionality.tsv"))[, qtl := "mqtl"]
he <- fread(file.path(PATHS$results, "coloc_eqtl", "hits_directionality.tsv"))[, qtl := "eqtl"]
hits <- rbind(hm, he, fill = TRUE)
snps <- unique(hits$top_snp)
count_step("Hits to audit", nrow(hits))

# GWAS
g <- fread(PATHS$gwas_raw, skip = "CHROM", select = c("CHROM", "ID", "POS", "A1", "A2", "BETA"))
g <- g[ID %in% snps]
setnames(g, c("CHR", "SNP", "POS", "gA1", "gA2", "gBETA"))
g[, `:=`(gA1 = toupper(gA1), gA2 = toupper(gA2))]
g[, raw_risk := fifelse(gBETA > 0, gA1, gA2)]
hits <- merge(hits, g, by.x = "top_snp", by.y = "SNP", all.x = TRUE)
count_step("Top SNPs found in the raw GWAS file", sum(!is.na(hits$gBETA)))

# mQTL
m <- fread(PATHS$mqtl_raw, select = c("SNP_Chr", "SNP_BP", "SNP_Allele", "ProbeID", "beta"))
m <- m[ProbeID %in% hits[qtl == "mqtl"]$feature]
m[, SNP_Allele := toupper(SNP_Allele)]
m <- m[nchar(SNP_Allele) == 1]   # SNVs only
hits <- merge(hits, m, by.x = c("feature", "CHR", "POS"), by.y = c("ProbeID", "SNP_Chr", "SNP_BP"), all.x = TRUE)

# eQTL
sub_e <- file.path(PATHS$processed, "eqtl_loci_subset.txt")
bed   <- file.path(PATHS$processed, "snp_positions_subset.bed")
if (!file.exists(sub_e) || !file.exists(bed)) stop("eQTL subset files from 04_qc_eqtl.R not found in ", PATHS$processed)
e <- fread(sub_e, select = c("gene_id", "variant_id", "slope"))[gene_id %in% hits[qtl == "eqtl"]$feature & variant_id %in% snps]
b <- fread(bed, header = FALSE, select = c(4, 7, 8), col.names = c("variant_id", "REF", "ALT"))[variant_id %in% snps]
b <- b[!duplicated(variant_id)]
e <- merge(e, b, by = "variant_id")
hits <- merge(hits, e, by.x = c("feature", "top_snp"), by.y = c("gene_id", "variant_id"), all.x = TRUE)

# comparison
hits[, raw_qtl_effect_of_risk := fifelse(qtl == "mqtl",
         fifelse(SNP_Allele == raw_risk, beta, fifelse(SNP_Allele == gA1 | SNP_Allele == gA2, -beta, NA_real_)),
         fifelse(ALT == raw_risk & (REF == gA1 | REF == gA2), slope,
                 fifelse(REF == raw_risk & (ALT == gA1 | ALT == gA2), -slope, NA_real_)))]
hits[, raw_direction := fifelse(is.na(raw_qtl_effect_of_risk), NA_character_,
         paste(fifelse(raw_qtl_effect_of_risk > 0, "increases", "decreases"),
               fifelse(qtl == "mqtl", "methylation", "expression")))]
hits[, risk_allele_match := risk_allele == raw_risk]
hits[, direction_match   := effect_of_risk_allele == raw_direction]
hits[, status := fifelse(is.na(raw_direction), "could not audit (raw row missing)",
                 fifelse(risk_allele_match & direction_match, "match", "MISMATCH"))]

out <- hits[, .(qtl, locus_id, feature, feature_label, top_snp,
                gwas_A1 = gA1, gwas_A2 = gA2, gwas_BETA_A1 = gBETA, raw_risk_allele = raw_risk,
                mqtl_SNP_Allele = SNP_Allele, mqtl_beta = beta, eqtl_REF = REF, eqtl_ALT = ALT, eqtl_slope_ALT = slope,
                raw_qtl_effect_of_risk, raw_direction,
                pipeline_risk_allele = risk_allele, pipeline_direction = effect_of_risk_allele, status)]
fwrite(out[order(qtl, status, locus_id)], file.path(PATHS$results, "direction_audit.tsv"), sep = "\t")
for (s in unique(out$status)) count_step(paste("Audit status:", s), sum(out$status == s))
options(width = 220)
print(out[status != "match"])
finish_log()
