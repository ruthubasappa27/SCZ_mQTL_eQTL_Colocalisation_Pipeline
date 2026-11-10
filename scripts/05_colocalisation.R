# 05_colocalisation.R
# Colocalisation of the schizophrenia GWAS with fetal brain mQTL or eQTL: one test per
# feature (CpG or gene) per locus.
#
# Each test runs coloc.abf and coloc.susie (SuSiE fine-mapping of each trait). The primary
# result is coloc.susie when both traits have >= 1 credible set, otherwise coloc.abf; the
# reason is recorded. For coloc.susie, the primary result is the credible-set pair with the
# highest PP4, and the same pair is followed in the p12 sensitivity analysis.
#
# SNPs whose allele frequency in the GWAS or QTL data differs from the LD panel by more than
# SETTINGS$max_freq_diff are removed before each test (QTL frequencies are available for the eQTL
# data only); the number removed is recorded.
#
# Per-locus results are cached with a fingerprint (MD5 of the input files, LD panel, PLINK
# binary, settings, result-producing code and package versions); a locus is recomputed
# whenever its fingerprint differs.
#
# Output: results/coloc_<qtl>/coloc_all_tests.tsv, coloc_snpwise_pp.tsv.gz, failed_loci.tsv
# Usage:  Rscript scripts/05_colocalisation.R mqtl|eqtl [n_loci]

source("config/config.R")
source("scripts/functions.R")
suppressMessages({ library(coloc); library(susieR) })

args  <- commandArgs(trailingOnly = TRUE)
qtl   <- if (length(args) >= 1) args[1] else "mqtl"
if (!qtl %in% c("mqtl", "eqtl")) stop("Usage: Rscript scripts/05_colocalisation.R mqtl|eqtl [n_loci]")
limit <- if (length(args) >= 2) as.integer(args[2]) else NA

start_log(paste0("05_colocalisation_", qtl))
outdir <- file.path(PATHS$results, paste0("coloc_", qtl))
perdir <- file.path(outdir, "per_locus")
dir.create(perdir, showWarnings = FALSE, recursive = TRUE)

gwas  <- fread(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"))
q     <- fread(file.path(PATHS$processed, paste0(qtl, "_harmonised.tsv.gz")))
loci  <- fread(file.path(PATHS$results, "loci", "loci_for_coloc.tsv"))
if (!is.na(limit)) loci <- head(loci, limit)

n_qtl <- SETTINGS[[paste0("n_", qtl)]]
s_cc  <- SETTINGS$gwas_case_fraction   # N is the effective (balanced) sample size
P1 <- SETTINGS$p1; P2 <- SETTINGS$p2; P12 <- SETTINGS$p12
PPcols <- paste0("PP.H", 0:4, ".abf")

check_susie_n_purity()
quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressWarnings(suppressMessages(expr))))
  out
}

# one colocalisation test (one locus x one feature)
run_test <- function(L, fi, gl, qf, ref) {
  base <- data.table(locus_id = L$locus_id, feature = fi$FEATURE,
                     feature_label = fi$FEATURE_LABEL, feature_chr = fi$FEATURE_CHR,
                     feature_pos = fi$FEATURE_POS)
  snps <- Reduce(intersect, list(gl$SNP, qf$SNP, rownames(ref$R)))
  fr   <- ref$freq_a1[snps]
  snps <- snps[is.finite(fr) & fr > 0 & fr < 1]
  freq_bad <- function(f) (abs(f - ref$freq_a1[snps]) > SETTINGS$max_freq_diff) %in% TRUE
  bad <- freq_bad(gl$FREQ_A1[match(snps, gl$SNP)])
  if ("FREQ_A1" %in% names(qf)) bad <- bad | freq_bad(qf$FREQ_A1[match(snps, qf$SNP)])
  n_freq <- sum(bad)
  snps <- snps[!bad]
  if (length(snps) < SETTINGS$min_overlap_snps)
    return(list(row = cbind(base, status = "not_tested", nsnps = length(snps), n_freq_mismatch = n_freq,
                            reason = sprintf("only %d SNPs shared by GWAS, QTL and LD panel", length(snps)))))

  g1  <- gl[match(snps, SNP)]
  q1  <- qf[match(snps, SNP)]
  maf <- unname(pmin(ref$freq_a1[snps], 1 - ref$freq_a1[snps]))
  d1  <- list(snp = snps, position = g1$POS, beta = g1$BETA, varbeta = g1$SE^2,
              type = "cc", s = s_cc, N = round(median(g1$N)))
  d2  <- list(snp = snps, position = q1$POS, beta = q1$BETA, varbeta = q1$SE^2,
              type = "quant", N = n_qtl, MAF = maf)

  abf <- quiet(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12))

  # SuSiE fine-mapping of each trait, then coloc.susie
  R <- ref$R[snps, snps]
  d1L <- c(d1, list(LD = R)); d2L <- c(d2, list(LD = R))
  set.seed(SETTINGS$seed)
  s1 <- tryCatch(quiet(runsusie(d1L, L = SETTINGS$susie_max_signals,
                                n_purity = SETTINGS$susie_n_purity)), error = function(e) e)
  s2 <- tryCatch(quiet(runsusie(d2L, L = SETTINGS$susie_max_signals,
                                n_purity = SETTINGS$susie_n_purity)), error = function(e) e)
  err <- c(if (inherits(s1, "error")) paste("GWAS:", conditionMessage(s1)),
           if (inherits(s2, "error")) paste("QTL:",  conditionMessage(s2)))
  ncs1 <- if (inherits(s1, "error")) NA_integer_ else length(s1$sets$cs)
  ncs2 <- if (inherits(s2, "error")) NA_integer_ else length(s2$sets$cs)

  cs_res <- NULL
  if (length(err) == 0 && ncs1 > 0 && ncs2 > 0) {
    cs_res <- tryCatch(quiet(coloc.susie(s1, s2, p1 = P1, p2 = P2, p12 = P12)), error = function(e) NULL)
    if (!is.null(cs_res) && (is.null(cs_res$summary) || nrow(cs_res$summary) == 0 ||
                             all(is.na(cs_res$summary$PP.H4.abf)))) cs_res <- NULL
  }

  if (!is.null(cs_res)) {
    method <- "coloc.susie"
    reason <- "SuSiE credible set(s) in both traits"
    k      <- which.max(cs_res$summary$PP.H4.abf)
    pair_cols <- if (all(c("idx1", "idx2") %in% names(cs_res$summary))) c("idx1", "idx2")
                 else if (all(c("hit1", "hit2") %in% names(cs_res$summary))) c("hit1", "hit2")
                 else stop("coloc.susie summary has no signal-pair identifiers (idx1/idx2 or hit1/hit2)")
    pair   <- as.data.table(cs_res$summary)[k, ..pair_cols]
    susie_pair <- if (all(c("hit1", "hit2") %in% names(cs_res$summary)))
                    paste0(cs_res$summary$hit1[k], " | ", cs_res$summary$hit2[k]) else
                    paste0("cs", pair[[1]], " | cs", pair[[2]])
    pp_col <- if (nrow(cs_res$summary) == 1) "SNP.PP.H4.abf" else paste0("SNP.PP.H4.row", k)
    snp_pp <- data.table(snp = cs_res$results$snp, pp = cs_res$results[[pp_col]])
    prim   <- unlist(cs_res$summary[k, ..PPcols])
  } else {
    method <- "coloc.abf"
    reason <- if (length(err)) paste("SuSiE error -", paste(err, collapse = "; "))
              else if (isTRUE(ncs1 > 0 && ncs2 > 0)) "coloc.susie returned no comparable credible-set pair"
              else sprintf("no SuSiE credible set (GWAS: %s, QTL: %s)", ncs1, ncs2)
    snp_pp <- data.table(snp = abf$results$snp, pp = abf$results$SNP.PP.H4)
    prim   <- abf$summary[PPcols]
    susie_pair <- NA_character_
  }
  cs <- credible_set(snp_pp$snp, snp_pp$pp, SETTINGS$credible_coverage)

  # PP4 of the primary method across p12 values
  sens <- sapply(SETTINGS$sensitivity_p12, function(p12) {
    tryCatch({
      if (method == "coloc.susie") {
        sm  <- as.data.table(quiet(coloc.susie(s1, s2, p1 = P1, p2 = P2, p12 = p12))$summary)
        hit <- sm[pair, on = pair_cols, nomatch = NULL]
        if (nrow(hit) != 1) stop("selected credible-set pair not found at p12 = ", p12)
        hit$PP.H4.abf
      } else quiet(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = p12))$summary[["PP.H4.abf"]]
    }, error = function(e) NA_real_)
  })
  names(sens) <- paste0("PP4_p12_", format(SETTINGS$sensitivity_p12, scientific = TRUE))

  # LD mismatch between summary statistics and panel (estimate_s_rss)
  s_g <- tryCatch(estimate_s_rss(d1$beta / sqrt(d1$varbeta), R, n = d1$N), error = function(e) NA_real_)
  s_q <- tryCatch(estimate_s_rss(d2$beta / sqrt(d2$varbeta), R, n = n_qtl), error = function(e) NA_real_)

  row <- cbind(base, data.table(
    status = "tested", reason = reason, nsnps = length(snps), n_freq_mismatch = n_freq,
    gwas_minp_overlap = min(g1$P), qtl_minp_overlap = min(q1$P),
    gwas_lead_in_overlap = L$lead_snp %in% snps,
    abf_PP.H0 = abf$summary[["PP.H0.abf"]], abf_PP.H1 = abf$summary[["PP.H1.abf"]],
    abf_PP.H2 = abf$summary[["PP.H2.abf"]], abf_PP.H3 = abf$summary[["PP.H3.abf"]],
    abf_PP.H4 = abf$summary[["PP.H4.abf"]],
    susie_ncs_gwas = ncs1, susie_ncs_qtl = ncs2,
    susie_PP.H4_max = if (!is.null(cs_res)) max(cs_res$summary$PP.H4.abf) else NA_real_,
    method = method, susie_pair = susie_pair,
    PP.H0 = prim[[1]], PP.H1 = prim[[2]], PP.H2 = prim[[3]], PP.H3 = prim[[4]], PP4 = prim[[5]],
    top_snp = cs$snp[1], top_snp_pp = cs$pp[1], cs95_size = nrow(cs),
    ld_mismatch_s_gwas = s_g, ld_mismatch_s_qtl = s_q),
    as.data.table(as.list(sens)))
  snp_pp[, `:=`(locus_id = L$locus_id, feature = fi$FEATURE, method = method)]
  list(row = row, snp_pp = snp_pp)
}

# one locus: all features with a QTL association in the window
run_locus <- function(L) {
  none <- function(why) list(rows = data.table(locus_id = L$locus_id, feature = NA_character_,
                                               status = "not_tested", reason = why))
  gl <- gwas[CHR == L$CHR & POS >= L$start & POS <= L$end]
  ql <- q[CHR == L$CHR & POS >= L$start & POS <= L$end]
  if (nrow(gl) == 0) return(none("no harmonised GWAS SNPs in locus"))
  if (nrow(ql) == 0) return(none("no QTL associations in locus"))

  feats <- ql[, .(minp = min(P), FEATURE_LABEL = FEATURE_LABEL[1],
                  FEATURE_CHR = FEATURE_CHR[1], FEATURE_POS = FEATURE_POS[1]), by = FEATURE]
  thr  <- SETTINGS$qtl_feature_min_p
  rows <- list(); pps <- list()
  weak <- feats[minp >= thr]
  if (nrow(weak) > 0)
    rows[[1]] <- weak[, .(locus_id = L$locus_id, feature = FEATURE, feature_label = FEATURE_LABEL,
                          feature_chr = FEATURE_CHR, feature_pos = FEATURE_POS, status = "not_tested",
                          reason = sprintf("no QTL SNP with p < %g in locus", thr))]
  test <- feats[minp < thr]
  if (nrow(test) == 0) return(list(rows = rbindlist(rows, fill = TRUE)))

  snps_all <- intersect(gl$SNP, ql[FEATURE %in% test$FEATURE]$SNP)
  if (length(snps_all) < SETTINGS$min_overlap_snps) {
    rows[[length(rows) + 1]] <- test[, .(locus_id = L$locus_id, feature = FEATURE,
                                         feature_label = FEATURE_LABEL, feature_chr = FEATURE_CHR,
                                         feature_pos = FEATURE_POS, status = "not_tested",
                                         nsnps = length(snps_all),
                                         reason = "fewer shared GWAS-QTL SNPs than minimum")]
    return(list(rows = rbindlist(rows, fill = TRUE)))
  }

  ref <- plink_ld(L$CHR, snps_all, tag = L$locus_id)
  for (i in seq_len(nrow(test))) {
    fi <- test[i]
    r <- tryCatch(run_test(L, fi, gl, ql[FEATURE == fi$FEATURE], ref),
                  error = function(e) list(row = data.table(locus_id = L$locus_id, feature = fi$FEATURE,
                                                            feature_label = fi$FEATURE_LABEL,
                                                            status = "error", reason = conditionMessage(e))))
    rows[[length(rows) + 1]] <- r$row
    pps[[length(pps) + 1]]   <- r$snp_pp
  }
  list(rows = rbindlist(rows, fill = TRUE), snp_pp = rbindlist(pps, fill = TRUE))
}

# run fingerprint
fp_files <- c(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"),
              file.path(PATHS$processed, paste0(qtl, "_harmonised.tsv.gz")),
              file.path(PATHS$results, "loci", "loci_for_coloc.tsv"),
              as.vector(outer(paste0(PATHS$ld_prefix, 1:22), c(".bed", ".bim", ".fam"), paste0)),
              PATHS$plink)
if (!all(file.exists(fp_files))) stop("Missing input file(s): ", paste(fp_files[!file.exists(fp_files)], collapse = ", "))
fingerprint <- md5_of_text(c(
  paste(basename(fp_files), unname(tools::md5sum(fp_files))),
  qtl, n_qtl,
  deparse(SETTINGS[c("gwas_case_fraction", "n_mqtl", "n_eqtl", "qtl_feature_min_p", "min_overlap_snps", "max_freq_diff",
                     "p1", "p2", "p12", "susie_max_signals", "susie_n_purity", "seed",
                     "credible_coverage", "sensitivity_p12")]),
  unlist(lapply(list(run_test, run_locus, credible_set, plink_ld, load_bim),
                function(fn) deparse(utils::removeSource(fn)))),
  R.version.string, paste("coloc", packageVersion("coloc")), paste("susieR", packageVersion("susieR")),
  paste("data.table", packageVersion("data.table"))))
logmsg("  Run fingerprint: ", fingerprint)
is_current <- function(f) file.exists(f) && identical(readRDS(f)$fingerprint, fingerprint)

# main loop
count_step("Loci in this run", nrow(loci))
failed <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]
  f <- file.path(perdir, paste0(L$locus_id, ".rds"))
  if (is_current(f)) next
  t0  <- Sys.time()
  res <- tryCatch(run_locus(L), error = function(e) e)
  if (inherits(res, "error")) {
    logmsg(sprintf("[%d/%d] %s  ERROR: %s", i, nrow(loci), L$locus_id, conditionMessage(res)))
    failed[[length(failed) + 1]] <- data.table(locus_id = L$locus_id, error = conditionMessage(res))
    next
  }
  res$fingerprint <- fingerprint
  saveRDS(res, f)
  tested <- res$rows[status == "tested"]
  logmsg(sprintf("[%d/%d] %s  tests: %d  max PP4: %s  (%.0fs)", i, nrow(loci), L$locus_id,
                 nrow(tested), if (nrow(tested)) sprintf("%.3f", max(tested$PP4)) else "-",
                 as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# combine results
done <- lapply(file.path(perdir, paste0(loci$locus_id, ".rds")),
               function(f) if (is_current(f)) readRDS(f))
tests  <- rbindlist(lapply(done, `[[`, "rows"), fill = TRUE)
snp_pp <- rbindlist(lapply(done, `[[`, "snp_pp"), fill = TRUE)
fwrite(tests,  file.path(outdir, "coloc_all_tests.tsv"), sep = "\t")
fwrite(snp_pp, file.path(outdir, "coloc_snpwise_pp.tsv.gz"), sep = "\t")
fwrite(if (length(failed)) rbindlist(failed) else data.table(locus_id = character(), error = character()),
       file.path(outdir, "failed_loci.tsv"), sep = "\t")

count_step("Loci finished", sum(!sapply(done, is.null)))
count_step("Loci that stopped with an error (see failed_loci.tsv)", length(failed))
count_step("Loci with >=1 colocalisation test", uniqueN(tests[status == "tested"]$locus_id))
count_step("Colocalisation tests run", sum(tests$status == "tested"))
count_step("  of which coloc.susie primary", sum(tests$method == "coloc.susie", na.rm = TRUE))
count_step("  of which coloc.abf primary", sum(tests$method == "coloc.abf", na.rm = TRUE))
count_step("Tests with an error (see coloc_all_tests.tsv)", sum(tests$status == "error"))
count_step(sprintf("Tests with PP4 >= %.1f", SETTINGS$pp4_threshold),
           sum(tests$PP4 >= SETTINGS$pp4_threshold, na.rm = TRUE))
finish_log()
