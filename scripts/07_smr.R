# 07_smr.R
# SMR and HEIDI for every colocalised feature, at the primary instrument threshold
# (p < 5e-8) and a sensitivity threshold (p < 1e-4). SMR p-values are Bonferroni-corrected
# over features with an instrument; HEIDI p > 0.05 is a pass, and NA (fewer than three SNPs)
# is reported as not computable.
#
# mQTL: SMR-format Hannon et al. data (FB_Brain_2; SNP IDs are chr:bp).
# eQTL: BESD built from the harmonised O'Brien et al. data. Gene TSS positions are converted
#       from GRCh38 to GRCh37 with the Ensembl REST assembly mapping; each mapping is checked
#       (one mapping, same chromosome, single base) and the script stops if any fails.
# GWAS: .ma file with freq = frequency of A1.
#
# Output: results/smr_<qtl>/smr_results.tsv (one row per feature and analysis)
#         results/smr_eqtl/eqtl_gene_tss_GRCh38_to_GRCh37.tsv
# Usage:  Rscript scripts/07_smr.R mqtl|eqtl

source("config/config.R")
source("scripts/functions.R")
suppressMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
qtl  <- if (length(args) >= 1) args[1] else "mqtl"
if (!qtl %in% c("mqtl", "eqtl")) stop("Usage: Rscript scripts/07_smr.R mqtl|eqtl")
start_log(paste0("07_smr_", qtl))
outdir <- file.path(PATHS$results, paste0("smr_", qtl))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$smr_ref, showWarnings = FALSE, recursive = TRUE)

# helpers
# delete all files with this output prefix
clear_prefix <- function(prefix) {
  old <- list.files(dirname(prefix), pattern = paste0("^", basename(prefix), "\\."), full.names = TRUE)
  if (length(old)) unlink(old)
}

# run SMR; stop with the end of the console output on a non-zero exit status
run_smr <- function(args, prefix) {
  log_file <- paste0(prefix, ".console.txt")
  status <- system2(PATHS$smr, args, stdout = log_file, stderr = log_file)
  if (!identical(as.integer(status), 0L)) {
    tail_txt <- if (file.exists(log_file)) paste(tail(readLines(log_file, warn = FALSE), 15), collapse = "\n") else ""
    stop("SMR failed (exit status ", status, "). Last lines of ", log_file, ":\n", tail_txt)
  }
  invisible(log_file)
}

norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

# input checks before writing SMR files
assert_assoc <- function(d, what, freq_col) {
  bad <- c(
    "missing or non-SNV allele"  = sum(!(d$A1 %in% c("A","C","G","T")) | !(d$A2 %in% c("A","C","G","T"))),
    "A1 equals A2"               = sum(d$A1 == d$A2, na.rm = TRUE),
    "beta not finite"            = sum(!is.finite(d$BETA)),
    "SE not finite or <= 0"      = sum(!is.finite(d$SE) | d$SE <= 0),
    "P not in (0, 1]"            = sum(!is.finite(d$P) | d$P <= 0 | d$P > 1),
    "frequency not in (0, 1)"    = sum(!is.finite(d[[freq_col]]) | d[[freq_col]] <= 0 | d[[freq_col]] >= 1),
    "chromosome not 1-22"        = sum(!(as.character(d$CHR) %in% as.character(1:22))))
  if (any(bad > 0))
    stop(what, " failed input checks: ", paste(sprintf("%s = %d rows", names(bad)[bad > 0], bad[bad > 0]), collapse = "; "))
  logmsg("  ", what, ": all ", nrow(d), " rows pass input checks (alleles, beta, SE, P, frequency, chromosome)")
}

# colocalised features
hits <- fread(file.path(PATHS$results, paste0("coloc_", qtl), "hits.tsv"))
if (nrow(hits) == 0) stop("No colocalised features to test.")
probes <- unique(hits$feature)
count_step("Colocalised features to test with SMR", length(probes))
multi <- hits[, .N, by = feature][N > 1]
count_step("Features colocalised in more than one locus (listed together in one row)", nrow(multi))
coloc_by_feature <- hits[order(-PP4), .(
  feature_label  = paste(unique(na.omit(feature_label)), collapse = ";"),
  n_coloc_hits   = .N,
  coloc_hits     = paste(sprintf("%s (PP4 %.3f)", locus_id, PP4), collapse = "; ")),
  by = .(probeID = feature)]

# QTL data in BESD format
if (qtl == "mqtl") {
  besd <- PATHS$besd_mqtl
  for (ext in c(".besd", ".epi", ".esi"))
    if (!file.exists(paste0(besd, ext))) stop("Missing mQTL BESD file: ", besd, ext)
  esi <- fread(paste0(besd, ".esi"), header = FALSE, nrows = 1000)
  use_chrbp <- !any(grepl("^rs", esi$V2))
} else {
  use_chrbp <- FALSE
  q <- fread(file.path(PATHS$processed, "eqtl_harmonised.tsv.gz"))[FEATURE %in% probes]
  missing_genes <- setdiff(probes, unique(q$FEATURE))
  assert_assoc(q, "eQTL associations for SMR", "FREQ_A1")
  if (anyDuplicated(q[, .(FEATURE, SNP)])) stop("Duplicate gene-SNP rows in the eQTL data for SMR")
  if (length(missing_genes)) stop("Colocalised genes absent from eqtl_harmonised.tsv.gz: ", paste(missing_genes, collapse = ", "))

  # GRCh38 TSS per gene (one chromosome and position per gene)
  tss <- q[, .(n_chr = uniqueN(CHR[!is.na(CHR)]), n_pos = uniqueN(FEATURE_POS[!is.na(FEATURE_POS)]),
               CHR = CHR[!is.na(CHR)][1], TSS38 = FEATURE_POS[!is.na(FEATURE_POS)][1]), by = FEATURE]
  bad <- tss[n_chr != 1 | n_pos != 1]
  if (nrow(bad)) stop("Genes without one consistent chromosome/TSS in the eQTL data: ", paste(bad$FEATURE, collapse = ", "))
  tss[, c("n_chr", "n_pos") := NULL]
  tss[, CHR := as.integer(norm_chr(CHR))]

  # GRCh38 -> GRCh37 with the Ensembl REST assembly mapping
  map_tss <- function(chr, pos) {
    chr <- norm_chr(chr); pos <- as.integer(pos)
    url <- sprintf("https://rest.ensembl.org/map/human/GRCh38/%s:%d..%d:1/GRCh37?content-type=application/json",
                   chr, pos, pos)
    body <- tempfile(fileext = ".json"); hdr <- tempfile(fileext = ".txt")
    on.exit(unlink(c(body, hdr)), add = TRUE)
    code <- NA_character_
    for (attempt in 1:5) {
      code <- tryCatch(system2("curl", c("-sS", "--max-time", "30", "-D", shQuote(hdr), "-o", shQuote(body),
                                         "-w", shQuote("%{http_code}"), shQuote(url)),
                               stdout = TRUE, stderr = FALSE), error = function(e) "000")
      code <- if (length(code)) tail(code, 1) else "000"
      if (code == "200") break
      if (code == "429") {                         # rate limited: wait for Retry-After
        h <- if (file.exists(hdr)) readLines(hdr, warn = FALSE) else character()
        ra <- suppressWarnings(as.numeric(trimws(sub("(?i)^retry-after:", "", grep("(?i)^retry-after:", h, value = TRUE, perl = TRUE), perl = TRUE))))
        Sys.sleep(if (length(ra) && is.finite(ra[1])) ra[1] + 1 else 2^attempt)
      } else if (code %in% c("000", "500", "502", "503", "504")) {
        Sys.sleep(2^attempt)                         # network or server error: retry
      } else break
    }
    out <- list(TSS38 = pos, TSS37 = NA_integer_, mapped_chr = NA_character_, mapped_strand = NA_integer_,
                http_status = code, map_status = NA_character_)
    if (code != "200") { out$map_status <- paste0("http_error_", code); return(out) }
    js <- tryCatch(fromJSON(body, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(js) || is.null(js$mappings)) { out$map_status <- "parse_error"; return(out) }
    m <- js$mappings
    if (length(m) == 0) { out$map_status <- "unmapped"; return(out) }
    if (length(m) > 1)  { out$map_status <- "split_mapping"; return(out) }
    o <- m[[1]]$original; mp <- m[[1]]$mapped
    out$mapped_chr <- norm_chr(mp$seq_region_name); out$mapped_strand <- as.integer(mp$strand)
    if (!identical(as.integer(o$start), pos))          { out$map_status <- "original_mismatch"; return(out) }
    if (!identical(out$mapped_chr, chr))                { out$map_status <- "other_chromosome"; return(out) }
    if (!identical(as.integer(mp$start), as.integer(mp$end))) { out$map_status <- "not_single_base"; return(out) }
    # strand of the mapped segment does not change a single-base coordinate; it is recorded
    if (!(out$mapped_strand %in% c(1L, -1L)))            { out$map_status <- "unknown_strand"; return(out) }
    if (!is.null(mp$assembly) && !identical(as.character(mp$assembly), "GRCh37")) { out$map_status <- "wrong_assembly"; return(out) }
    out$TSS37 <- as.integer(mp$start)
    out$map_status <- if (out$mapped_strand == 1L) "ok" else "ok_reverse_strand_segment"
    out
  }
  mapped <- rbindlist(lapply(seq_len(nrow(tss)), function(i) {
    Sys.sleep(0.1)
    as.data.table(map_tss(tss$CHR[i], tss$TSS38[i]))[, FEATURE := tss$FEATURE[i]]
  }))
  tss <- merge(tss[, .(FEATURE, CHR)], mapped, by = "FEATURE")
  fwrite(tss, file.path(outdir, "eqtl_gene_tss_GRCh38_to_GRCh37.tsv"), sep = "\t")
  for (st in unique(tss$map_status)) count_step(paste("eQTL gene TSS mapping:", st), sum(tss$map_status == st))
  bad <- tss[!map_status %in% c("ok", "ok_reverse_strand_segment")]
  if (nrow(bad))
    stop("GRCh38 -> GRCh37 TSS mapping not valid for ", nrow(bad), " gene(s): ",
         paste(sprintf("%s (%s)", bad$FEATURE, bad$map_status), collapse = ", "),
         "\nSee ", file.path(outdir, "eqtl_gene_tss_GRCh38_to_GRCh37.tsv"))

  # build the BESD file
  edir <- file.path(outdir, "esd")
  unlink(edir, recursive = TRUE); dir.create(edir)
  for (f in tss$FEATURE)
    fwrite(q[FEATURE == f, .(Chr = CHR, SNP, Bp = POS, A1, A2, Freq = FREQ_A1, Beta = BETA, se = SE, p = P)],
           file.path(edir, paste0(f, ".esd")), sep = "\t")
  flist <- tss[, .(Chr = CHR, ProbeID = FEATURE, GeneticDistance = 0, ProbeBp = TSS37,
                   Gene = FEATURE, Orientation = "NA", PathOfEsd = file.path(edir, paste0(FEATURE, ".esd")))]
  fwrite(flist, file.path(outdir, "eqtl.flist"), sep = "\t")
  besd <- file.path(outdir, "OBrien_coloc_genes")
  clear_prefix(besd)
  run_smr(c("--eqtl-flist", file.path(outdir, "eqtl.flist"), "--make-besd", "--out", besd), besd)
  for (ext in c(".besd", ".epi", ".esi"))
    if (!file.exists(paste0(besd, ext))) stop("Building the eQTL BESD did not produce ", besd, ext,
                                              ": see ", besd, ".console.txt")
}
logmsg("  BESD SNP IDs are ", if (use_chrbp) "chr:bp" else "rsIDs")
epi <- fread(paste0(besd, ".epi"), header = FALSE,
             col.names = c("CHR", "probe", "cm", "bp", "gene", "strand"))
epi <- epi[probe %in% probes]
epi[, CHR := suppressWarnings(as.integer(norm_chr(CHR)))]
if (any(is.na(epi$CHR) | !epi$CHR %in% 1:22)) stop("BESD probe chromosomes are not all 1-22: check ", besd, ".epi")
count_step("Colocalised features present in the BESD file", nrow(epi))

# GWAS summary statistics in .ma format
g <- fread(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"))
if (use_chrbp) g[, SNP := paste0(CHR, ":", POS)]
g <- g[!SNP %in% SNP[duplicated(SNP)]]
assert_assoc(g, "GWAS summary statistics for SMR", "FREQ_A1")
ma <- file.path(outdir, "gwas_for_smr.ma")
fwrite(g[, .(SNP, A1, A2, freq = FREQ_A1, b = BETA, se = SE, p = P, n = round(N))], ma, sep = "\t")

# LD reference with SNP IDs matching the BESD file
# A chr:bp copy of the panel is reused only if its stamp records the same source files.
ld_ref <- function(chr) {
  src <- paste0(PATHS$ld_prefix, chr)
  src_files <- paste0(src, c(".bed", ".bim", ".fam"))
  if (!all(file.exists(src_files))) stop("Missing LD panel file(s) for chr", chr, ": ",
                                         paste(src_files[!file.exists(src_files)], collapse = ", "))
  if (!use_chrbp) return(src)
  out   <- file.path(PATHS$smr_ref, paste0("EUR_chr", chr, "_chrbp"))
  stamp <- paste0(out, ".source_md5.txt")
  src_md5 <- paste(basename(src_files), unname(tools::md5sum(src_files)))
  trio_ok <- all(file.exists(paste0(out, c(".bed", ".bim", ".fam"))))
  if (trio_ok && file.exists(stamp) && identical(readLines(stamp, warn = FALSE), src_md5)) return(out)

  clear_prefix(out)
  bim <- fread(paste0(src, ".bim"), header = FALSE)
  id  <- paste0(bim$V1, ":", bim$V4)
  dup <- bim$V2[id %in% id[duplicated(id)]]          # positions shared by >1 SNP
  args <- c("--bfile", src, "--keep-allele-order", "--make-bed", "--out", out)
  if (length(dup)) { writeLines(dup, paste0(out, ".dups")); args <- c(args, "--exclude", paste0(out, ".dups")) }
  run_plink(args)
  if (!all(file.exists(paste0(out, c(".bed", ".bim", ".fam"))))) stop("PLINK did not produce the full trio for ", out)
  b2 <- fread(paste0(out, ".bim"), header = FALSE)
  b2[, V2 := paste0(V1, ":", V4)]
  if (anyDuplicated(b2$V2)) stop("Duplicate chr:bp IDs remain in ", out, ".bim")
  fwrite(b2, paste0(out, ".bim"), sep = "\t", col.names = FALSE)
  writeLines(src_md5, stamp)
  out
}

# SMR per chromosome and threshold
thresholds <- c(primary = SETTINGS$smr_peqtl_main, sensitivity = SETTINGS$smr_peqtl_sensitivity)
res <- list()
for (chr in sort(unique(epi$CHR))) {
  pf <- file.path(outdir, paste0("probes_chr", chr, ".txt"))
  writeLines(epi[CHR == chr]$probe, pf)
  ref <- ld_ref(chr)
  for (an in names(thresholds)) {
    out <- file.path(outdir, sprintf("smr_%s_chr%s", an, chr))
    clear_prefix(out)
    run_smr(c("--bfile", ref, "--gwas-summary", ma,
              "--beqtl-summary", besd, "--extract-probe", pf,
              "--peqtl-smr", thresholds[[an]], "--thread-num", SETTINGS$smr_threads,
              "--out", out), out)
    if (!file.exists(paste0(out, ".smr")))
      stop("SMR exited without error but wrote no ", out, ".smr; see ", out, ".console.txt")
    r <- fread(paste0(out, ".smr"))
    if (nrow(r)) res[[length(res) + 1]] <- cbind(analysis = an, r)
  }
}
res <- rbindlist(res, fill = TRUE)
for (col in c("analysis", "probeID", "topSNP", "p_SMR", "p_HEIDI", "nsnp_HEIDI"))
  if (!col %in% names(res)) res[, (col) := if (col %in% c("analysis", "probeID", "topSNP")) character() else numeric()]
if (anyDuplicated(res[, .(analysis, probeID)])) stop("SMR returned more than one row for a probe in one analysis")

# one row per probe and analysis, including probes without an instrument
grid <- CJ(analysis = names(thresholds), probeID = probes)
res  <- merge(grid, res, by = c("analysis", "probeID"), all.x = TRUE)
res[, peqtl := thresholds[analysis]]
res[, n_tested := sum(!is.na(p_SMR)), by = analysis]
res[, `:=`(smr_bonferroni_sig = !is.na(p_SMR) & p_SMR < 0.05 / n_tested,
           heidi_pass = fifelse(is.na(p_HEIDI), NA, p_HEIDI > 0.05),
           note = fifelse(!is.na(p_SMR), fifelse(is.na(p_HEIDI), "HEIDI not computable (inconclusive)", ""),
                  fifelse(probeID %in% epi$probe, "no instrument at this p threshold", "probe not in BESD")))]
res <- merge(res, coloc_by_feature, by = "probeID", all.x = TRUE)
if (nrow(res) != length(probes) * length(thresholds)) stop("Unexpected number of rows in the SMR results table")
fwrite(res, file.path(outdir, "smr_results.tsv"), sep = "\t")

for (an in names(thresholds)) {
  r <- res[analysis == an]
  count_step(sprintf("SMR %s (p_eQTL < %g): probes with instrument", an, thresholds[[an]]), sum(!is.na(r$p_SMR)))
  count_step(sprintf("SMR %s: Bonferroni-significant", an), sum(r$smr_bonferroni_sig))
  count_step(sprintf("SMR %s: significant AND HEIDI p > 0.05", an), sum(r$smr_bonferroni_sig & r$heidi_pass %in% TRUE))
  count_step(sprintf("SMR %s: significant, HEIDI not computable", an), sum(r$smr_bonferroni_sig & is.na(r$heidi_pass)))
}
options(width = 200)
print(res[, .(analysis, probeID, feature_label, n_coloc_hits, topSNP, p_SMR, p_HEIDI, nsnp_HEIDI, note)])
finish_log()
