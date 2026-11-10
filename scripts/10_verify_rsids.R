# 10_verify_rsids.R
# Checks every reported SNP against Ensembl GRCh37 (REST, POST /variation/homo_sapiens):
# the colocalisation top SNP of every colocalised feature and the SMR top SNP of every
# SMR-tested feature. SMR chr:bp IDs (mQTL) are converted to rsIDs with the 1000G EUR panel
# (same position and alleles). For each rsID: found, current ID if merged, same chromosome
# and position as the panel, and both panel alleles in Ensembl's allele string.
# Status: OK, OK_merged, not_found, not_in_panel, position_mismatch, allele_mismatch.
#
# Output: results/rsid_check/rsid_check.tsv, smr_top_snp_rsid.tsv
# Usage:  Rscript scripts/10_verify_rsids.R

source("config/config.R")
source("scripts/functions.R")
suppressMessages(library(jsonlite))
start_log("10_verify_rsids")
outdir <- file.path(PATHS$results, "rsid_check")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bim <- load_bim()

# SNP IDs
used <- list()
for (q in c("mqtl", "eqtl")) {
  h <- fread(file.path(PATHS$results, paste0("coloc_", q), "hits.tsv"))
  used[[length(used) + 1]] <- h[, .(qtl = q, feature, role = "coloc_top_snp", snp_reported = top_snp)]
  s <- fread(file.path(PATHS$results, paste0("smr_", q), "smr_results.tsv"))[analysis == "primary" & !is.na(topSNP) & topSNP != ""]
  if (nrow(s)) used[[length(used) + 1]] <- s[, .(qtl = q, feature = probeID, role = "smr_top_snp",
                                                snp_reported = topSNP, smr_chr = topSNP_chr, smr_bp = topSNP_bp,
                                                smr_A1 = A1, smr_A2 = A2)]
}
used <- rbindlist(used, fill = TRUE)

# chr:bp IDs -> rsIDs (same position and allele pair in the panel)
cb <- used[!grepl("^rs[0-9]+$", snp_reported)]
if (nrow(cb)) {
  m <- merge(cb[, .(snp_reported, CHR = as.integer(smr_chr), POS = as.integer(smr_bp), smr_A1, smr_A2)],
             bim[, .(CHR, POS, rsid = SNP, A1, A2)], by = c("CHR", "POS"))
  m <- unique(m[(A1 == smr_A1 & A2 == smr_A2) | (A1 == smr_A2 & A2 == smr_A1), .(snp_reported, rsid)])
  amb <- m[, .N, by = snp_reported][N > 1]$snp_reported
  m <- m[!snp_reported %in% amb]
  used[m, rsid := i.rsid, on = "snp_reported"]
  count_step("SMR chr:bp IDs", uniqueN(cb$snp_reported))
  count_step("SMR chr:bp IDs converted to an rsID", uniqueN(m$snp_reported))
  count_step("SMR chr:bp IDs with >1 matching panel rsID (left unconverted)", length(amb))
}
used[grepl("^rs[0-9]+$", snp_reported), rsid := snp_reported]
fwrite(used[role == "smr_top_snp", .(qtl, feature, snp_reported, rsid)],
       file.path(outdir, "smr_top_snp_rsid.tsv"), sep = "\t")

ids <- sort(unique(na.omit(used$rsid)))
count_step("Unique rsIDs to check", length(ids))

# Ensembl GRCh37 lookup, 200 IDs per request
post_ensembl <- function(batch) {
  body <- tempfile(fileext = ".json"); out <- tempfile(fileext = ".json"); hdr <- tempfile()
  on.exit(unlink(c(body, out, hdr)))
  writeLines(toJSON(list(ids = batch)), body)
  for (attempt in 1:6) {
    code <- system2("curl", c("-sS", "--max-time", "60", "-X", "POST",
                              "-H", shQuote("Content-Type: application/json"),
                              "-H", shQuote("Accept: application/json"),
                              "-D", shQuote(hdr), "-o", shQuote(out), "-w", shQuote("%{http_code}"),
                              "--data", paste0("@", shQuote(body)),
                              shQuote("https://grch37.rest.ensembl.org/variation/homo_sapiens")),
                    stdout = TRUE, stderr = FALSE)
    code <- tail(code, 1)
    if (identical(code, "200")) return(fromJSON(out, simplifyVector = FALSE))
    wait <- if (code == "429") {
      ra <- grep("^retry-after:", tolower(readLines(hdr, warn = FALSE)), value = TRUE)
      if (length(ra)) as.numeric(sub(".*:\\s*", "", ra[1])) + 1 else 5
    } else 2^attempt
    logmsg(sprintf("  Ensembl HTTP %s, retrying in %.0f s", code, wait)); Sys.sleep(wait)
  }
  stop("Ensembl REST request failed after 6 attempts (last HTTP status ", code, ").")
}
resp <- list()
for (b in split(ids, ceiling(seq_along(ids) / 200))) resp <- c(resp, post_ensembl(b))

# comparison with the panel
comp <- function(a) c(A = "T", C = "G", G = "C", T = "A")[a]
pb <- bim[SNP %in% ids]
chk <- rbindlist(lapply(ids, function(id) {
  p <- pb[SNP == id]
  if (nrow(p) != 1) return(data.table(rsid = id, status = "not_in_panel"))
  r <- resp[[id]]
  if (is.null(r)) return(data.table(rsid = id, status = "not_found"))
  mp <- Filter(function(x) identical(x$assembly_name, "GRCh37") || is.null(x$assembly_name), r$mappings)
  pos_ok <- any(vapply(mp, function(x) isTRUE(as.character(x$seq_region_name) == as.character(p$CHR)) &&
                                      isTRUE(as.integer(x$start) == p$POS), logical(1)))
  al <- unique(unlist(lapply(mp, function(x) strsplit(x$allele_string, "/")[[1]])))
  al_ok <- all(c(p$A1, p$A2) %in% al) || all(comp(c(p$A1, p$A2)) %in% al)
  current <- if (is.null(r$name)) id else r$name
  status <- if (!pos_ok) "position_mismatch" else if (!al_ok) "allele_mismatch" else
            if (current != id) "OK_merged" else "OK"
  data.table(rsid = id, status = status, current_rsid = current,
             panel_chr = p$CHR, panel_pos = p$POS, panel_alleles = paste(p$A1, p$A2, sep = "/"),
             ensembl_grch37 = paste(vapply(mp, function(x) paste0(x$seq_region_name, ":", x$start), ""), collapse = ";"),
             ensembl_alleles = paste(al, collapse = "/"))
}), fill = TRUE)

out <- merge(used, chk, by = "rsid", all.x = TRUE)
out[is.na(rsid), status := "no_rsid (chr:bp not converted)"]
fwrite(out[order(status, qtl, role, feature)], file.path(outdir, "rsid_check.tsv"), sep = "\t")
for (st in sort(unique(chk$status))) count_step(sprintf("rsIDs with status %s", st), sum(chk$status == st))
bad <- chk[!status %in% c("OK", "OK_merged")]
if (nrow(bad)) { logmsg("  rsIDs needing attention:"); print(bad) } else logmsg("  All rsIDs verified against Ensembl GRCh37.")
finish_log()
