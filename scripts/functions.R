# functions.R
# Helper functions shared by the pipeline scripts.

suppressMessages(library(data.table))

# logging
.log_env <- new.env()

start_log <- function(name) {
  dir.create(PATHS$logs, showWarnings = FALSE, recursive = TRUE)
  .log_env$file   <- file.path(PATHS$logs, paste0(name, ".log"))
  .log_env$name   <- name
  .log_env$counts <- data.table(script = character(), step = character(), n = numeric())
  cat("", file = .log_env$file)
  logmsg("== ", name, " started ", format(Sys.time(), "%Y-%m-%d %H:%M"), " ==")
}

logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  if (!is.null(.log_env$file)) cat(msg, "\n", file = .log_env$file, append = TRUE)
}

# record a count (logged and written to results/qc_counts.tsv)
count_step <- function(step, n) {
  logmsg(sprintf("  %-60s %s", step, format(n, big.mark = ",")))
  .log_env$counts <- rbind(.log_env$counts,
                           data.table(script = .log_env$name, step = step, n = n))
}

finish_log <- function() {
  dir.create(PATHS$results, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(PATHS$results, "qc_counts.tsv")
  old <- if (file.exists(f)) fread(f)[script != .log_env$name] else NULL
  fwrite(rbind(old, .log_env$counts), f, sep = "\t")
  writeLines(capture.output(sessionInfo()),
             file.path(PATHS$logs, paste0(.log_env$name, "_sessionInfo.txt")))
  logmsg("== ", .log_env$name, " finished ", format(Sys.time(), "%Y-%m-%d %H:%M"), " ==")
}

need_cols <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0)
    stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
         "\nColumns found: ", paste(names(dt), collapse = ", "))
}

# alleles
is_snv <- function(a) a %in% c("A", "C", "G", "T")
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
  (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# reference panel
load_bim <- function(chrs = 1:22) {
  bim <- rbindlist(lapply(chrs, function(chr) {
    fread(paste0(PATHS$ld_prefix, chr, ".bim"), header = FALSE,
          col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"),
          colClasses = list(character = c(2, 5, 6)))
  }))
  bim[, CHR := as.integer(CHR)]
  bim <- bim[is_snv(A1) & is_snv(A2)]
  # remove rsIDs that occur more than once
  bim <- bim[!SNP %in% SNP[duplicated(SNP)]]
  bim
}

# PLINK
run_plink <- function(args) {
  status <- system2(PATHS$plink, c(args, "--silent"), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("PLINK failed. Arguments: ", paste(args, collapse = " "))
  invisible(TRUE)
}

# signed LD (r, coded on panel A1) and A1 frequency for a set of SNPs
plink_ld <- function(chr, snps, tag) {
  dir.create(PATHS$tmp, showWarnings = FALSE, recursive = TRUE)
  pre <- file.path(PATHS$tmp, tag)
  writeLines(unique(snps), paste0(pre, ".snps"))
  run_plink(c("--bfile", paste0(PATHS$ld_prefix, chr),
              "--extract", paste0(pre, ".snps"),
              "--keep-allele-order", "--r", "square", "--write-snplist",
              "--out", pre))
  run_plink(c("--bfile", paste0(PATHS$ld_prefix, chr),
              "--extract", paste0(pre, ".snps"),
              "--keep-allele-order", "--freq", "--out", pre))
  ids <- readLines(paste0(pre, ".snplist"))
  R <- as.matrix(fread(paste0(pre, ".ld"), header = FALSE))
  dimnames(R) <- list(ids, ids)
  # drop SNPs with undefined LD (monomorphic in the panel)
  repeat {
    bad <- rownames(R)[apply(R, 1, function(r) any(!is.finite(r)))]
    if (length(bad) == 0) break
    worst <- names(which.max(colSums(!is.finite(R[, bad, drop = FALSE]))))
    R <- R[rownames(R) != worst, colnames(R) != worst, drop = FALSE]
  }
  frq <- fread(paste0(pre, ".frq"))
  freq <- setNames(frq$MAF, frq$SNP)   # A1 frequency (--keep-allele-order)
  unlink(list.files(PATHS$tmp, pattern = paste0("^", tag, "\\."), full.names = TRUE))
  list(R = R, freq_a1 = freq)
}

# PDF device that keeps Unicode symbols (quartz on macOS, cairo elsewhere, base pdf as fallback)
pdf_device <- function() {
  if (Sys.info()[["sysname"]] == "Darwin" && capabilities("aqua"))
    return(function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height))
  ok <- tryCatch({ f <- tempfile(fileext = ".pdf"); grDevices::cairo_pdf(f); grDevices::dev.off(); unlink(f); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (ok) grDevices::cairo_pdf else "pdf"
}

# reproducibility
md5_of_text <- function(x) {
  tf <- tempfile(); on.exit(unlink(tf)); writeLines(x, tf); unname(tools::md5sum(tf))
}

# Checks that susieR computes credible-set purity over all SNPs when n_purity = Inf.
# Test data: 150 SNPs with z = 5 and r = 0.9, except r = 0.85 between SNPs 1 and 2; the
# credible set (> 100 SNPs) includes both, so its exact purity is 0.85.
check_susie_n_purity <- function() {
  p <- 150; R <- matrix(0.9, p, p); diag(R) <- 1; R[1, 2] <- R[2, 1] <- 0.85
  fit <- tryCatch(suppressWarnings(suppressMessages(
           susieR::susie_rss(z = rep(5, p), R = R, n = 10000, L = 1,
                             n_purity = SETTINGS$susie_n_purity))),
         error = function(e) e)
  if (inherits(fit, "error"))
    stop("Installed susieR does not accept n_purity: ", conditionMessage(fit))
  cs <- if (length(fit$sets$cs)) fit$sets$cs[[1]] else integer(0)
  if (length(cs) <= 100 || !all(1:2 %in% cs))
    stop("susieR purity self-test could not build the test credible set (size ", length(cs),
         ").")
  exact <- min(abs(R[cs, cs])); got <- fit$sets$purity$min.abs.corr[1]
  if (!isTRUE(all.equal(got, exact, tolerance = 1e-8)))
    stop(sprintf("susieR purity self-test failed: reported %.6f, exact %.6f.",
                 got, exact))
  logmsg(sprintf("  susieR purity self-test passed (n_purity = %s; exact purity %.2f reproduced)",
                 SETTINGS$susie_n_purity, exact))
  invisible(TRUE)
}

# gene names
# Current HGNC symbol (org.Hs.eg.db). keytype "SYMBOL": a current symbol is kept, an outdated
# symbol is replaced through a unique alias match; keytype "ENSEMBL": unique Ensembl -> symbol.
# Inputs without a unique match are returned unchanged (symbols) or as the Ensembl ID.
current_symbol <- function(x, keytype = c("SYMBOL", "ENSEMBL")) {
  keytype <- match.arg(keytype)
  suppressMessages(requireNamespace("org.Hs.eg.db"))
  db <- org.Hs.eg.db::org.Hs.eg.db
  out <- x
  ok <- !is.na(x) & x != ""
  if (!any(ok)) return(out)
  one <- function(keys, kt, col) {
    r <- tryCatch(suppressMessages(AnnotationDbi::mapIds(db, keys, col, kt, multiVals = "list")),
                  error = function(e) setNames(vector("list", length(keys)), keys))
    vapply(keys, function(k) { v <- unique(na.omit(unlist(r[[k]]))); if (length(v) == 1) v else NA_character_ }, "")
  }
  if (keytype == "ENSEMBL") {
    sym <- one(sub("\\..*$", "", x[ok]), "ENSEMBL", "SYMBOL")
    out[ok] <- ifelse(is.na(sym), sub("\\..*$", "", x[ok]), sym)
  } else {
    valid <- x[ok] %in% AnnotationDbi::keys(db, "SYMBOL")
    alias <- one(x[ok][!valid], "ALIAS", "SYMBOL")
    tmp <- x[ok]; tmp[!valid] <- ifelse(is.na(alias), x[ok][!valid], alias)
    out[ok] <- tmp
  }
  out
}

# credible sets
credible_set <- function(snp, pp, coverage = 0.95) {
  o <- order(pp, decreasing = TRUE)
  dt <- data.table(snp = snp[o], pp = pp[o])
  dt[, cum_pp := cumsum(pp) / sum(pp)]
  w <- which(dt$cum_pp >= coverage)[1]   # include the SNP that crosses the threshold
  dt[seq_len(w)]
}
