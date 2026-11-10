# 06_results.R
# Summarises the colocalisation results of one QTL dataset: locus overview, colocalised
# features (PP4 >= threshold), 95% credible sets with LD to the top SNP and to the GWAS lead
# SNPs, direction of effect of the schizophrenia risk allele, and stability across p12.
#
# Output: results/coloc_<qtl>/locus_overview.tsv, hits.tsv, hits_credible_sets.tsv,
#         hits_directionality.tsv
# Usage:  Rscript scripts/06_results.R mqtl|eqtl

source("config/config.R")
source("scripts/functions.R")
args <- commandArgs(trailingOnly = TRUE)
qtl  <- if (length(args) >= 1) args[1] else "mqtl"
start_log(paste0("06_results_", qtl))
outdir <- file.path(PATHS$results, paste0("coloc_", qtl))
thr    <- SETTINGS$pp4_threshold
trait  <- if (qtl == "mqtl") "methylation" else "expression"

tests  <- fread(file.path(outdir, "coloc_all_tests.tsv"))
snp_pp <- fread(file.path(outdir, "coloc_snpwise_pp.tsv.gz"))
loci   <- fread(file.path(PATHS$results, "loci", "loci_all.tsv"))
ff     <- file.path(outdir, "failed_loci.tsv")
failed <- if (file.exists(ff) && isTRUE(file.size(ff) > 0)) fread(ff) else data.table(locus_id = character())

# locus overview
best_of <- function(x, by) as.character(if (length(by) && any(!is.na(by))) x[which.max(by)] else NA)
lo <- tests[, .(
  n_features_tested   = sum(status == "tested"),
  n_features_untested = sum(status == "not_tested"),
  n_errors            = sum(status == "error"),
  max_PP4      = if (any(status == "tested")) max(PP4[status == "tested"]) else NA_real_,
  best_feature = best_of(feature[status == "tested"], PP4[status == "tested"]),
  why_not_tested = paste(unique(reason[status != "tested"]), collapse = " | ")
), by = locus_id]
lo <- merge(loci[, .(locus_id, CHR, start, end, lead_snp, lead_p, n_lead_snps, status)],
            lo, by = "locus_id", all.x = TRUE)
lo[locus_id %in% failed$locus_id, status := "pipeline_error"]
lo[status == "analysed" & is.na(n_features_tested), status := "not_run_yet"]
lo[status == "analysed", status := fifelse(n_features_tested > 0, "usable", "not_usable")]
lo[, colocalised := !is.na(max_PP4) & max_PP4 >= thr]
fwrite(lo, file.path(outdir, "locus_overview.tsv"), sep = "\t")
for (s in unique(lo$status)) count_step(paste("Loci with status:", s), sum(lo$status == s))
count_step(sprintf("Loci with >=1 feature PP4 >= %.1f", thr), sum(lo$colocalised))

# colocalised features
hits <- tests[status == "tested" & PP4 >= thr][order(-PP4)]
hits[, feature_label := sapply(strsplit(fifelse(is.na(feature_label), "", feature_label), ";"),
                              function(x) paste(unique(x[nzchar(x)]), collapse = ";"))]
count_step(sprintf("Colocalised %s features (PP4 >= %.1f)", qtl, thr), nrow(hits))
count_step("  in distinct loci", uniqueN(hits$locus_id))
# descriptive: genome-wide significant GWAS SNP among the shared SNPs
hits[, gwas_gws_in_overlap := gwas_minp_overlap < SETTINGS$gwas_p_sig]
count_step("  of which overlap contains a GWS (p<5e-8) GWAS SNP", sum(hits$gwas_gws_in_overlap))
count_step("  ... in distinct loci", uniqueN(hits[gwas_gws_in_overlap == TRUE]$locus_id))
if (nrow(hits) == 0) { fwrite(hits, file.path(outdir, "hits.tsv"), sep = "\t"); finish_log(); quit(save = "no") }

sens_cols <- grep("^PP4_p12_", names(hits), value = TRUE)
hits[, robust_to_p12 := apply(.SD, 1, function(x) all(x[1:2] >= thr, na.rm = FALSE)), .SDcols = sens_cols]
count_step("  of which PP4 >= threshold also at p12 = 1e-6 and 5e-6", sum(hits$robust_to_p12, na.rm = TRUE))

# credible sets and LD
gwas <- fread(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"))
q    <- fread(file.path(PATHS$processed, paste0(qtl, "_harmonised.tsv.gz")))
cs_all <- list(); dir_all <- list()
for (i in seq_len(nrow(hits))) {
  h  <- hits[i]
  L  <- loci[locus_id == h$locus_id]
  pp <- snp_pp[locus_id == h$locus_id & feature == h$feature]
  cs <- credible_set(pp$snp, pp$pp, SETTINGS$credible_coverage)
  leads <- strsplit(L$merged_leads, ",")[[1]]          # all lead SNPs in the locus
  ld <- tryCatch(plink_ld(L$CHR, unique(c(cs$snp, leads)), tag = paste0("cs_", i)),
                 error = function(e) NULL)
  r2 <- function(a, b) if (!is.null(ld) && a %in% rownames(ld$R) && b %in% rownames(ld$R)) ld$R[a, b]^2 else NA_real_
  cs[, `:=`(locus_id = h$locus_id, feature = h$feature, feature_label = h$feature_label,
            r2_with_top_snp  = sapply(snp, r2, b = cs$snp[1]),
            r2_with_gwas_lead = sapply(snp, function(x) suppressWarnings(max(sapply(leads, r2, a = x), na.rm = TRUE))),
            nearest_gwas_lead = sapply(snp, function(x) { v <- sapply(leads, r2, a = x)
                                         if (all(is.na(v))) NA_character_ else leads[which.max(v)] }))]
  cs[!is.finite(r2_with_gwas_lead), r2_with_gwas_lead := NA_real_]
  cs_all[[i]] <- cs

  # direction of effect at the top SNP
  top <- cs$snp[1]
  g1 <- gwas[SNP == top]; q1 <- q[FEATURE == h$feature & SNP == top]
  if (nrow(g1) == 1 && nrow(q1) == 1) {
    risk <- if (g1$BETA > 0) g1$A1 else g1$A2
    qtl_beta_risk <- q1$BETA * sign(g1$BETA)     # QTL effect of the risk allele
    dir_all[[i]] <- data.table(locus_id = h$locus_id, feature = h$feature, feature_label = h$feature_label,
      top_snp = top, risk_allele = risk, other_allele = if (risk == g1$A1) g1$A2 else g1$A1,
      gwas_OR_risk = exp(abs(g1$BETA)), gwas_logOR_risk = abs(g1$BETA), gwas_p = g1$P,
      qtl_beta_risk = qtl_beta_risk, qtl_p = q1$P,
      effect_of_risk_allele = paste(ifelse(qtl_beta_risk > 0, "increases", "decreases"), trait))
  }
}
cs_all  <- rbindlist(cs_all, fill = TRUE)
dir_all <- rbindlist(dir_all, fill = TRUE)
if (nrow(dir_all) > 0) hits <- merge(hits, dir_all[, .(locus_id, feature, risk_allele, qtl_beta_risk, effect_of_risk_allele)],
              by = c("locus_id", "feature"), all.x = TRUE)[order(-PP4)]

fwrite(hits,    file.path(outdir, "hits.tsv"), sep = "\t")
fwrite(cs_all,  file.path(outdir, "hits_credible_sets.tsv"), sep = "\t")
fwrite(dir_all, file.path(outdir, "hits_directionality.tsv"), sep = "\t")

# r2 between the top SNP and the nearest GWAS lead SNP of the locus
lead_r2 <- cs_all[, .(r2_top_vs_gwas_lead = r2_with_gwas_lead[1], gwas_lead = nearest_gwas_lead[1],
                      cs_max_r2_gwas_lead = suppressWarnings(max(r2_with_gwas_lead, na.rm = TRUE))),
                  by = .(locus_id, feature)]
lead_r2[!is.finite(cs_max_r2_gwas_lead), cs_max_r2_gwas_lead := NA_real_]
hits <- merge(hits, lead_r2, by = c("locus_id", "feature"), all.x = TRUE)[order(-PP4)]
fwrite(hits, file.path(outdir, "hits.tsv"), sep = "\t")
options(width = 250)
logmsg("\nColocalised signals (PP4 >= threshold). Descriptive columns: gwas_minp = smallest GWAS p among shared SNPs;")
logmsg("r2_lead = LD between top colocalised SNP and nearest GWAS lead SNP of the locus; robust_to_p12 = PP4 >= threshold at p12 1e-6 and 5e-6")
print(hits[, .(locus_id, feature, gene = feature_label, method, nsnps, PP4 = round(PP4, 3),
               gwas_minp = signif(gwas_minp_overlap, 2), gwas_lead, r2_lead = round(r2_top_vs_gwas_lead, 2),
               cs95_size, robust_to_p12)])
finish_log()
