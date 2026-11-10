# 04_qc_eqtl.R
# Quality control of the O'Brien et al. (2018) fetal brain eQTL data and alignment to the
# 1000 Genomes EUR reference panel. Only SNPs within the analysed loci are used.
#
# Input:  all_eqtls_gene.txt (slope per ALT allele) and snp_positions.bed (GRCh38 position,
#         rsID, REF, ALT, ALT frequency), from the O'Brien et al. data release
# Output: data/processed/eqtl_harmonised.tsv.gz
#         (FEATURE = Ensembl gene ID; FEATURE_POS = TSS, GRCh38; BETA = effect of panel allele A1)
# Usage:  Rscript scripts/04_qc_eqtl.R
#
# SNPs are matched to the GRCh37 panel by rsID, and REF/ALT must match the panel alleles.

source("config/config.R")
source("scripts/functions.R")
start_log("04_qc_eqtl")
dir.create(PATHS$tmp, showWarnings = FALSE, recursive = TRUE)

snp_bed <- PATHS$eqtl_snp_pos
if (!file.exists(snp_bed) && file.exists(paste0(snp_bed, ".gz"))) snp_bed <- paste0(snp_bed, ".gz")
if (!file.exists(snp_bed)) stop("data/snp_positions.bed not found (O'Brien figshare).")

# harmonised GWAS SNPs inside the analysed loci
loci <- fread(file.path(PATHS$results, "loci", "loci_for_coloc.tsv"))
gwas <- fread(file.path(PATHS$processed, "gwas_harmonised.tsv.gz"), select = c("SNP", "CHR", "POS"))
setkey(gwas, CHR, POS)
keep <- unique(rbindlist(lapply(seq_len(nrow(loci)), function(i)
  gwas[CHR == loci$CHR[i] & POS >= loci$start[i] & POS <= loci$end[i], .(SNP)]))$SNP)
count_step("eQTL: GWAS SNPs inside analysed loci (filter list)", length(keep))
keep_file <- file.path(PATHS$tmp, "eqtl_keep_snps.txt")
writeLines(keep, keep_file)
rm(gwas); invisible(gc())

# subset the eQTL file to those SNPs
subset_file <- file.path(PATHS$processed, "eqtl_loci_subset.txt")
if (!file.exists(subset_file)) {
  logmsg("  Subsetting all_eqtls_gene.txt to SNPs in the analysed loci")
  cmd <- sprintf("LC_ALL=C awk 'NR==FNR{k[$1]=1; next} FNR==1 || ($2 in k)' '%s' '%s' > '%s.part' && mv '%s.part' '%s'",
                 keep_file, PATHS$eqtl_raw, subset_file, subset_file, subset_file)
  if (system(cmd) != 0) stop("awk filtering of the eQTL file failed")
}
e <- fread(subset_file)
need_cols(e, c("gene_id", "variant_id", "tss_distance", "maf", "pval_nominal", "slope", "slope_se"), "eQTL file")
count_step("eQTL: gene-SNP pairs inside analysed loci", nrow(e))
count_step("eQTL: genes", uniqueN(e$gene_id))

# REF/ALT alleles from snp_positions.bed
bed_subset <- file.path(PATHS$processed, "snp_positions_subset.bed")
reader <- if (grepl("\\.gz$", snp_bed)) sprintf("gunzip -c '%s'", snp_bed) else sprintf("cat '%s'", snp_bed)
cmd <- sprintf("%s | LC_ALL=C awk 'NR==FNR{k[$1]=1; next} ($4 in k)' '%s' - > '%s'",
               reader, keep_file, bed_subset)
if (system(cmd) != 0) stop("filtering snp_positions.bed failed")
s <- fread(bed_subset, header = FALSE, select = c(1, 3, 4, 7, 8, 9),
           col.names = c("chr38", "pos38", "SNP", "REF", "ALT", "INFO"))
s <- s[SNP %in% unique(e$variant_id)]
s[, AF_ALT := as.numeric(sub(".*(^|;)AF=([^;]+).*", "\\2", INFO))][, INFO := NULL]
s <- s[!SNP %in% SNP[duplicated(SNP)]]
e <- merge(e, s, by.x = "variant_id", by.y = "SNP")
count_step("eQTL: pairs with REF/ALT from snp_positions.bed", nrow(e))

# QC
e <- e[is.finite(slope) & is.finite(slope_se) & slope_se > 0 &
       is.finite(pval_nominal) & pval_nominal > 0 & pval_nominal <= 1]
count_step("eQTL: valid slope/SE/p", nrow(e))
e <- e[is_snv(REF) & is_snv(ALT)];                   count_step("eQTL: biallelic SNVs", nrow(e))

# alignment to the reference panel (rsID, then alleles)
bim <- load_bim()
e <- merge(e, bim[, .(SNP, CHR, POS, A1, A2)], by.x = "variant_id", by.y = "SNP")
count_step("eQTL: rsID present in 1000G EUR panel", nrow(e))
e[, align := fifelse(ALT == A1 & REF == A2, "same",
             fifelse(ALT == A2 & REF == A1, "swapped", "mismatch"))]
count_step("eQTL: REF/ALT do not match panel (removed)", sum(e$align == "mismatch"))
e <- e[align != "mismatch"]
e[, BETA    := fifelse(align == "same", slope, -slope)]       # effect of panel A1
e[, FREQ_A1 := fifelse(align == "same", AF_ALT, 1 - AF_ALT)]  # A1 frequency in O'Brien samples
count_step("eQTL: ALT = panel A2 (slope sign flipped)", sum(e$align == "swapped"))

dup <- duplicated(e[, .(gene_id, variant_id)])
count_step("eQTL: duplicate gene-SNP records removed", sum(dup))
e <- e[!dup]
count_step("eQTL: final gene-SNP pairs", nrow(e))
count_step("eQTL: final genes", uniqueN(e$gene_id))
count_step("eQTL: final SNPs", uniqueN(e$variant_id))

fwrite(e[, .(FEATURE = gene_id, FEATURE_LABEL = gene_id,
             FEATURE_CHR = CHR, FEATURE_POS = pos38 - tss_distance,   # TSS, GRCh38
             SNP = variant_id, CHR, POS, A1, A2, EA_original = ALT,
             BETA, SE = slope_se, P = pval_nominal, FREQ_A1)],
       file.path(PATHS$processed, "eqtl_harmonised.tsv.gz"), sep = "\t")
finish_log()
