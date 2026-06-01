# CADD Score and Variant Consequence Annotation
# For all credible set SNPs across 6 colocalised loci
# Author: Ruthushree P Basappa
# Date: May 2026


library(httr)
library(jsonlite)
library(data.table)
library(here)

here::i_am("05_functional_validation/cadd_annotation.R")

dir.create(here("results", "functional_validation"), showWarnings = FALSE, recursive = TRUE)

gene_snps <- list(
  OPRD1 = c("rs61783570", "rs482692", "rs369247", "rs204054"),
  ERBB4 = c("rs7607363", "rs6735626", "rs7355673", "rs10182996",
            "rs6710653", "rs2128321"),
  ZNF389 = c("rs13217472", "rs9257140", "rs3131341", "rs3118359",
             "rs3132385", "rs3132389"),
  PIK3C2B = c("rs4633059", "rs4546666", "rs13270146",
              "rs59880252", "rs13257396", "rs145454152"),
  WDR55 = c("rs2563335", "rs2530240", "rs801189", "rs1089305",
            "rs801185", "rs2563306", "rs813897", "rs2563304",
            "rs3776129", "rs2563302", "rs801168", "rs801167",
            "rs809635", "rs2530233", "rs2530232", "rs6849",
            "rs702395", "rs702394", "rs801183", "rs801182",
            "rs801180", "rs801179", "rs2563283", "rs2563336",
            "rs59690131"),
  DYNC1I2 = c("rs4667693", "rs10184866", "rs10201419",
              "rs6745143", "rs2678155", "rs7597387")
)

get_cadd <- function(snp) {
  url <- paste0(
    "https://myvariant.info/v1/variant/", snp,
    "?assembly=hg19&fields=cadd.phred,cadd.consequence,cadd.gene.genename"
  )
  tryCatch({
    response <- GET(url)
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    cadd_scores <- data$cadd$phred
    consequences <- data$cadd$consequence
    genes <- data$cadd$gene$genename
    cadd_val <- cadd_scores[!is.na(cadd_scores)][1]
    cons_val <- consequences[!is.na(consequences)][1]
    gene_val <- genes[!is.na(genes)][1]

    data.frame(
      SNP = snp,
      CADD_score = ifelse(is.null(cadd_val) || length(cadd_val) == 0, NA, cadd_val),
      Consequence = ifelse(is.null(cons_val) || length(cons_val) == 0, NA, cons_val),
      Gene = ifelse(is.null(gene_val) || length(gene_val) == 0, NA, gene_val)
    )
  }, error = function(e) {
    data.frame(SNP = snp, CADD_score = NA, Consequence = NA, Gene = NA)
  })
}

all_results <- list()
for (gene in names(gene_snps)) {
  cat("Processing", gene, "...\n")
  snp_list <- gene_snps[[gene]]
  gene_results <- rbindlist(lapply(snp_list, get_cadd))
  gene_results$Gene_candidate <- gene
  all_results[[gene]] <- gene_results
  Sys.sleep(0.5)
}

final_results <- rbindlist(all_results)
print(final_results)

fwrite(
  final_results,
  here("results", "functional_validation", "credible_set_CADD_scores_complete.csv")
)

cat("Saved to results/functional_validation/credible_set_CADD_scores_complete.csv\n")

