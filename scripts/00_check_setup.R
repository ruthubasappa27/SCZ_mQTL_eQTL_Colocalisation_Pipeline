# 00_check_setup.R
# Checks that the required R packages, programs, reference panel and input files are
# present, and prints the first rows of each input file.
#
# Usage: Rscript scripts/00_check_setup.R 2>&1 | tee logs/00_check_setup.txt

source("config/config.R")
ok <- TRUE
say <- function(status, msg) cat(sprintf("[%s] %s\n", status, msg))

cat("\n1. R packages\n")
pkgs <- c("data.table", "here", "R.utils", "coloc", "susieR", "jsonlite", "ggplot2", "gtable",
          "flextable", "officer", "openxlsx", "GEOquery", "org.Hs.eg.db", "AnnotationDbi")
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) say("OK", paste(p, packageVersion(p)))
  else { say("MISSING", p); ok <- FALSE }
}

cat("\n2. Programs\n")
for (tool in c("plink", "smr")) {
  f <- PATHS[[tool]]
  if (!file.exists(f)) { say("MISSING", f); ok <- FALSE; next }
  out <- tryCatch(suppressWarnings(system2(f, "--help", stdout = TRUE, stderr = TRUE)),
                  error = function(e) conditionMessage(e))
  first <- grep("PLINK|SMR|version", out, value = TRUE, ignore.case = TRUE)
  say("FOUND", paste0(tool, ": ", if (length(first)) first[1] else paste(head(out, 2), collapse = " ")))
}

cat("\n3. LD reference panel\n")
need <- as.vector(outer(paste0(PATHS$ld_prefix, 1:22), c(".bed", ".bim", ".fam"), paste0))
miss <- need[!file.exists(need)]
if (length(miss) == 0) {
  say("OK", "all 66 files (chr1-22 .bed/.bim/.fam)")
  cat("   first lines of chr22 .bim:\n"); print(data.table::fread(paste0(PATHS$ld_prefix, "22.bim"), nrows = 3))
} else {
  say("MISSING", paste(length(miss), "of 66 files, e.g.", paste(basename(head(miss, 3)), collapse = ", ")))
  cat("   files in references/ld_reference:\n")
  print(head(list.files(dirname(PATHS$ld_prefix)), 10))
  ok <- FALSE
}

cat("\n4. Data files\n")
cat("   files in data/:\n"); print(list.files("data"))
peek <- function(f, ...) {
  if (!file.exists(f)) {
    if (file.exists(paste0(f, ".gz"))) say("ZIPPED", paste0("run: gunzip ", f, ".gz"))
    else say("MISSING", f)
    ok <<- FALSE
    return(invisible(FALSE))
  }
  say("OK", sprintf("%s (%.0f MB)", basename(f), file.size(f) / 1e6))
  x <- data.table::fread(f, nrows = 3, ...)
  cat("   columns:", paste(names(x), collapse = " | "), "\n")
  print(x)
  invisible(TRUE)
}

cat("\n-- GWAS --\n")
if (isTRUE(peek(PATHS$gwas_raw, skip = "CHROM"))) {
  hdr <- grep("^##", readLines(PATHS$gwas_raw, n = 200), value = TRUE)
  cat("   header lines:\n"); cat(head(hdr, 40), sep = "\n")
}

cat("\n-- mQTL --\n")
if (isTRUE(peek(PATHS$mqtl_raw))) {
  m <- data.table::fread(PATHS$mqtl_raw)
  cat("   rows:", nrow(m), "\n")
  if ("p.value" %in% names(m))
    cat("   p == 0:", sum(m$p.value == 0, na.rm = TRUE), "| smallest p:", min(m$p.value, na.rm = TRUE), "\n")
  if ("ProbeID" %in% names(m)) cat("   unique CpGs:", data.table::uniqueN(m$ProbeID), "\n")
  if ("SNP_Allele" %in% names(m)) { cat("   SNP_Allele values:\n"); print(head(sort(table(m$SNP_Allele), decreasing = TRUE), 8)) }
}

cat("\n-- eQTL --\n")
peek(PATHS$eqtl_raw)
bed <- PATHS$eqtl_snp_pos
if (file.exists(bed) || file.exists(paste0(bed, ".gz"))) say("OK", basename(bed)) else { say("MISSING", bed); ok <- FALSE }

cat("\n-- CoDEx exports (needed for step 13 only) --\n")
for (f in file.path(PATHS$codex_dir, c("codex_mqtl.csv", "codex_eqtl.csv")))
  say(if (file.exists(f)) "OK" else "MISSING", f)

cat("\n5. Hannon BESD (SMR format)\n")
besd <- paste0(PATHS$besd_mqtl, c(".besd", ".epi", ".esi"))
if (all(file.exists(besd))) {
  say("OK", "FB_Brain_2 .besd/.epi/.esi")
  print(data.table::fread(besd[3], nrows = 3, header = FALSE))
  if (file.exists(PATHS$smr)) {
    dir.create(PATHS$tmp, showWarnings = FALSE)
    q <- file.path(PATHS$tmp, "besd_query")
    system2(PATHS$smr, c("--beqtl-summary", PATHS$besd_mqtl, "--query", "1",
                         "--probe", "cg16329650", "--out", q), stdout = FALSE, stderr = FALSE)
    if (file.exists(paste0(q, ".txt"))) {
      qq <- data.table::fread(paste0(q, ".txt"))
      pcol <- grep("^p", names(qq), value = TRUE)[1]
      # A BESD file with full cis results also stores non-significant SNPs
      cat("   probe cg16329650: SNPs stored =", nrow(qq),
          "| SNPs with p > 1e-5 =", sum(qq[[pcol]] > 1e-5), "\n")
    } else say("NOTE", "SMR query of the BESD file produced no output")
  }
} else { say("MISSING", "FB_Brain_2 files (needed for SMR only)") }

cat("\nRESULT\n")
cat(if (ok) "All required items found.\n" else "Some items are missing (see above).\n")
