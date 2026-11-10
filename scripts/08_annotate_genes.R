# 08_annotate_genes.R
# Gene annotation of every colocalised feature (labels only; nothing is filtered).
#   mQTL: Illumina gene annotation provided with the Hannon et al. data; nearest gene and
#         nearest protein-coding gene by TSS (UCSC hg19 refGene); all genes with a TSS
#         within 50 kb of the CpG.
#   eQTL: Ensembl GRCh37 symbol, biotype and position of each Ensembl gene ID.
#   gene_symbol: current HGNC symbol (org.Hs.eg.db) of each label or Ensembl ID, used for
#         reporting in all later steps; Ensembl ID where no HGNC symbol exists.
#
# Output: results/coloc_<qtl>/hits_gene_annotation.tsv
# Usage:  Rscript scripts/08_annotate_genes.R (downloads refGene and queries Ensembl REST)

source("config/config.R")
source("scripts/functions.R")
suppressMessages(library(jsonlite))
start_log("08_annotate_genes")
dir.create(PATHS$annotation_dir, showWarnings = FALSE, recursive = TRUE)

# reference gene table (hg19 refGene)
rg_file <- file.path(PATHS$annotation_dir, "refGene_hg19.txt.gz")
if (!file.exists(rg_file)) {
  url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
  if (system2("curl", c("-sSfL", "-o", rg_file, url)) != 0) stop("Download of refGene failed: ", url)
}
rg <- fread(rg_file, header = FALSE, select = c(2, 3, 4, 5, 6, 13),
            col.names = c("acc", "chrom", "strand", "txStart", "txEnd", "gene"))
rg <- rg[grepl("^chr[0-9]+$", chrom)]
rg[, CHR := as.integer(sub("chr", "", chrom))]
rg[, TSS := fifelse(strand == "+", txStart, txEnd)]
rg[, coding := grepl("^NM_", acc)]
rg <- unique(rg[, .(gene, CHR, TSS, coding)])
count_step("refGene hg19 transcripts start sites (unique gene x TSS)", nrow(rg))

annotate_cpg <- function(chr, pos) {
  g <- rg[CHR == chr]
  g[, d := abs(TSS - pos)]
  near   <- g[which.min(d)]
  nearc  <- g[coding == TRUE][which.min(d)]
  within <- unique(g[d <= 50000][order(d)]$gene)
  data.table(nearest_gene = near$gene, nearest_gene_dist_bp = near$d,
             nearest_coding = if (nrow(nearc)) nearc$gene else NA_character_,
             nearest_coding_dist_bp = if (nrow(nearc)) nearc$d else NA_real_,
             genes_within_50kb = paste(within, collapse = ";"))
}

# mQTL
m <- fread(file.path(PATHS$results, "coloc_mqtl", "hits.tsv"))
m[, illumina_label := sapply(strsplit(fifelse(is.na(feature_label), "", feature_label), ";"),
                             function(x) paste(unique(x[nzchar(x)]), collapse = ";"))]
ann_m <- rbindlist(lapply(seq_len(nrow(m)), function(i)
  cbind(m[i, .(locus_id, feature, feature_chr, feature_pos, illumina_label, PP4)],
        annotate_cpg(m$feature_chr[i], m$feature_pos[i]))))
ann_m[, illumina_agrees_nearest := mapply(function(a, b) b %in% strsplit(a, ";")[[1]],
                                          illumina_label, nearest_gene)]
# current HGNC symbols of the Illumina label genes (used for reporting)
lab_genes <- unique(unlist(strsplit(ann_m$illumina_label, ";")))
lab_genes <- lab_genes[nzchar(lab_genes)]
lab_map <- setNames(current_symbol(lab_genes, "SYMBOL"), lab_genes)
ann_m[, gene_symbol := vapply(strsplit(illumina_label, ";"), function(g)
  paste(unique(lab_map[g[nzchar(g)]]), collapse = ";"), "")]
count_step("mQTL label genes renamed to a current HGNC symbol", sum(lab_map != names(lab_map)))
fwrite(ann_m, file.path(PATHS$results, "coloc_mqtl", "hits_gene_annotation.tsv"), sep = "\t")
count_step("mQTL hits annotated", nrow(ann_m))
count_step("mQTL hits with blank Illumina label", sum(ann_m$illumina_label == ""))
count_step("mQTL hits where Illumina label includes the nearest gene", sum(ann_m$illumina_agrees_nearest))

# eQTL: Ensembl GRCh37 lookup
e <- fread(file.path(PATHS$results, "coloc_eqtl", "hits.tsv"))
ids <- unique(e$feature)
body <- toJSON(list(ids = ids), auto_unbox = FALSE)
tmpj <- file.path(PATHS$tmp, "ensembl_lookup.json")
dir.create(PATHS$tmp, showWarnings = FALSE)
status <- system2("curl", c("-sS", "-X", "POST", "https://grch37.rest.ensembl.org/lookup/id",
                            "-H", shQuote("Content-Type: application/json"),
                            "-H", shQuote("Accept: application/json"),
                            "-d", shQuote(body), "-o", tmpj))
lk <- if (status == 0 && file.exists(tmpj)) fromJSON(tmpj) else list()
ann_e <- rbindlist(lapply(ids, function(id) {
  x <- lk[[id]]
  if (is.null(x)) return(data.table(feature = id, symbol = NA_character_, biotype = NA_character_,
                                    ens_chr = NA_character_, ens_start = NA_real_, ens_end = NA_real_,
                                    note = "not found in Ensembl GRCh37"))
  data.table(feature = id, symbol = x$display_name %||% NA_character_, biotype = x$biotype %||% NA_character_,
             ens_chr = x$seq_region_name %||% NA_character_, ens_start = x$start %||% NA_real_,
             ens_end = x$end %||% NA_real_, note = "")
}))
ann_e <- merge(e[, .(locus_id, feature, PP4)], ann_e, by = "feature", all.x = TRUE)[order(-PP4)]
ann_e[, gene_symbol := current_symbol(feature, "ENSEMBL")]
# no unique Ensembl -> symbol match: use the Ensembl GRCh37 name if it is (or is an alias of) a
# current HGNC symbol; otherwise keep the Ensembl ID
fb <- ann_e[grepl("^ENSG", gene_symbol) & !is.na(symbol)]
if (nrow(fb)) {
  cand <- current_symbol(fb$symbol, "SYMBOL")
  ok <- cand %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
  ann_e[fb$feature[ok], gene_symbol := cand[ok], on = "feature"]
}
fwrite(ann_e, file.path(PATHS$results, "coloc_eqtl", "hits_gene_annotation.tsv"), sep = "\t")
count_step("eQTL hits annotated", nrow(ann_e))
count_step("eQTL hits not found in Ensembl GRCh37", sum(is.na(ann_e$symbol)))
count_step("eQTL hits without an HGNC symbol (Ensembl ID used)", sum(grepl("^ENSG", ann_e$gene_symbol)))

options(width = 220)
print(ann_m[, .(feature, illumina_label, gene_symbol, nearest_gene, nearest_gene_dist_bp, nearest_coding)])
print(ann_e[, .(feature, symbol, gene_symbol, biotype, ens_chr, ens_start, PP4 = round(PP4, 3))])
finish_log()
