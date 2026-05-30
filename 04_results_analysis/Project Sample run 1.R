Project: Smaple run 1 
Author: Ruthu 

#installing the packages 

install.packages("coloc")

install.packages("data.table")

install.packages("dplyr")

install.packages("ggplot2")

install.packages("ieugwasr") 

library(coloc)

library(data.table)

library(dplyr)

library(ggplot2)

library(ieugwasr)

# Download Trubetskoy 2022 schizophrenia GWAS

gwas <- associations(variants = "all", id = "ieu-b-5099") #REMEMBER TO REQUST ASSESS 

> gwas <- fread("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv")

#Read the file 
head(gwas)

colnames(gwas)

nrow(gwas)

#Dowmload Hannan et.al Data 

mqtl <- fread("All_Imputed_BonfSignificant_mQTLs.csv")

#Read the file 

head(mqtl)


#Claculate SE for the same 

#Calculate Z score from p-value

mqtl <- mqtl %>%
  mutate(
    Z = qnorm(p.value / 2, lower.tail = FALSE),
    SE = abs(beta) / Z
  )

head(mqtl %>% select(beta, p.value, Z, SE))

summary(mqtl$SE)

#Download eQTL data

eqtl <- fread("DER-08b_hg19_eQTL.bonferroni.txt")

head(eqtl)

colnames(eqtl)

eqtl <- eqtl %>%
  mutate(SNP_chr = gsub("chr", "", SNP_chr))

#Calculate SE 

eqtl <- eqtl %>%
  mutate(
    Z = qnorm(nominal_pval / 2, lower.tail = FALSE),
    SE = abs(regression_slope) / Z
  )
#Renaming the Coloums to match the standard 

eqtl <- eqtl %>% rename(
  CHR = SNP_chr,
  POS = SNP_start,
  BETA_eqtl = regression_slope,
  P_eqtl = nominal_pval
)

# Convert CHR to integer to match GWAS
eqtl <- eqtl %>%
  mutate(CHR = as.integer(CHR))

#Veiw the file 

head(eqtl %>% select(CHR, POS, BETA_eqtl, SE, P_eqtl))

#Sanity check 

cat("GWAS SNPs:", nrow(gwas), "\n")

cat("mQTL SNPs:", nrow(mqtl), "\n")

cat("eQTL SNPs:", nrow(eqtl), "\n")

#QC specific to GWAS 

# Step 1 — Remove low quality imputed SNPs

gwas <- gwas %>% filter(IMPINFO >= 0.8)

# Step 2 — Calculate MAF

gwas <- gwas %>%
  mutate(MAF = pmin(FCAS, 1 - FCAS))

# Step 3 — Remove rare variants

gwas <- gwas %>% filter(MAF >= 0.01)

# Step 4 — Keep only genome-wide significant SNPs

gwas_sig <- gwas %>% filter(P_gwas < 5e-8)

cat("SNPs after QC:", nrow(gwas), "\n")
cat("Genome-wide significant SNPs:", nrow(gwas_sig), "\n")

colnames(gwas)

# Rename columns to standard format

gwas <- gwas %>% rename(
  CHR = CHROM,
  SNP = ID,
  EA = A1,
  OA = A2,
  BETA_gwas = BETA,
  SE_gwas = SE,
  P_gwas = PVAL
)
#Rerunning 

gwas <- gwas %>% filter(IMPINFO >= 0.8)

gwas <- gwas %>% mutate(MAF = pmin(FCAS, 1 - FCAS))

gwas <- gwas %>% filter(MAF >= 0.01)

gwas_sig <- gwas %>% filter(P_gwas < 5e-8)

cat("SNPs after QC:", nrow(gwas), "\n")

cat("Genome-wide significant SNPs:", nrow(gwas_sig), "\n")

# Install susieR package

install.packages("susieR")

library(susieR)

# Install COLOC 

install.packages("coloc")

library(coloc)

packageVersion("coloc")


define_loci_fast <- function(gwas_sig, window = 1000000) {
  df <- gwas_sig %>% arrange(P_gwas)
  lead_snps <- data.frame()
  
  while(nrow(df) > 0) {
    top <- df[1, ]
    lead_snps <- rbind(lead_snps, top)
    df <- df %>%
      filter(!(CHR == top$CHR & 
                 abs(POS - top$POS) <= window/2))
  }
  return(lead_snps)
}

lead_snps <- define_loci_fast(gwas_sig)

cat("Number of independent loci:", nrow(lead_snps), "\n")

i <- 1
locus_chr <- lead_snps$CHR[i]
locus_start <- lead_snps$POS[i] - 500000
locus_end <- lead_snps$POS[i] + 500000

cat("Testing locus", i, "\n")
cat("Chr:", locus_chr, "Start:", locus_start, "End:", locus_end, "\n")

gwas_locus <- gwas %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)

mqtl_locus <- mqtl %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)

cat("GWAS SNPs in locus:", nrow(gwas_locus), "\n")
cat("mQTL SNPs in locus:", nrow(mqtl_locus), "\n")


  
colnames(gwas)

colnames(lead_snps)

i <- 1
locus_chr <- lead_snps$CHR[i]
locus_start <- lead_snps$POS[i] - 500000
locus_end <- lead_snps$POS[i] + 500000

cat("Chr:", locus_chr, "Start:", locus_start, "End:", locus_end, "\n")

gwas_locus <- gwas %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)

mqtl_locus <- mqtl %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)

cat("GWAS SNPs in locus:", nrow(gwas_locus), "\n")
cat("mQTL SNPs in locus:", nrow(mqtl_locus), "\n")

locus_chr <- lead_snps$CHR[1]
print(locus_chr)

locus_start <- lead_snps$POS[1] - 500000
print(locus_start)

locus_end <- lead_snps$POS[1] + 500000
print(locus_end)

gwas_locus <- gwas %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
cat("GWAS SNPs in locus:", nrow(gwas_locus), "\n")  

mqtl_locus <- mqtl %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
cat("mQTL SNPs in locus:", nrow(mqtl_locus), "\n")

colnames(mqtl)

mqtl <- mqtl %>% rename(
  CHR = SNP_Chr,
  POS = SNP_BP,
  EA_mqtl = SNP_Allele,
  BETA_mqtl = beta,
  P_mqtl = p.value,
  SE_mqtl = SE
)

colnames(mqtl)

mqtl_locus <- mqtl %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
cat("mQTL SNPs in locus:", nrow(mqtl_locus), "\n")

colnames(eqtl)
eqtl <- eqtl %>% rename(SE_eqtl = SE)

eqtl_locus <- eqtl %>%
  filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)

cat("eQTL SNPs in locus:", nrow(eqtl_locus), "\n")

#Merging the dataset 

#GWAS and mQTL 

merged_mqtl <- inner_join(gwas_locus, mqtl_locus, by = c("CHR", "POS"))
cat("GWAS + mQTL overlapping SNPs:", nrow(merged_mqtl), "\n")

#GWAS and eQTL

merged_eqtl <- inner_join(gwas_locus, eqtl_locus, by = c("CHR", "POS"))
cat("GWAS + eQTL overlapping SNPs:", nrow(merged_eqtl), "\n")

#Running COLOC for mQtl 

result_mqtl <- coloc.abf(
  dataset1 = list(
    beta = merged_mqtl$BETA_gwas,
    varbeta = merged_mqtl$SE_gwas^2,
    type = "cc",
    s = merged_mqtl$NCAS[1] / (merged_mqtl$NCAS[1] + merged_mqtl$NCON[1]),
    N = merged_mqtl$NCAS[1] + merged_mqtl$NCON[1],
    snp = merged_mqtl$SNP
  ),
  dataset2 = list(
    beta = merged_mqtl$BETA_mqtl,
    varbeta = merged_mqtl$SE_mqtl^2,
    type = "quant",
    N = 166,
    snp = merged_mqtl$SNP
  )
)

print(result_mqtl$summary)


#Running the loop 

results_mqtl <- list()
results_eqtl <- list()

for(i in 1:nrow(lead_snps)) {
  
  cat("Processing locus", i, "of", nrow(lead_snps), "\n")
  locus_chr <- lead_snps$CHR[i]
  locus_start <- lead_snps$POS[i] - 500000
  locus_end <- lead_snps$POS[i] + 500000
  lead_snp <- lead_snps$SNP[i]
  gwas_locus <- gwas %>%
    filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
  
  mqtl_locus <- mqtl %>%
    filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
  
  eqtl_locus <- eqtl %>%
    filter(CHR == locus_chr, POS >= locus_start, POS <= locus_end)
  if(nrow(gwas_locus) >= 100 & nrow(mqtl_locus) >= 100) {
    merged_mqtl <- inner_join(gwas_locus, mqtl_locus, by = c("CHR", "POS"))
    merged_mqtl <- merged_mqtl %>%
      group_by(SNP) %>%
      slice_min(P_mqtl, n = 1) %>%
      ungroup()
    
    if(nrow(merged_mqtl) >= 100) {
      
      tryCatch({
        
        result <- coloc.abf(
          dataset1 = list(
            beta = merged_mqtl$BETA_gwas,
            varbeta = merged_mqtl$SE_gwas^2,
            type = "cc",
            s = merged_mqtl$NCAS[1] / (merged_mqtl$NCAS[1] + merged_mqtl$NCON[1]),
            N = merged_mqtl$NCAS[1] + merged_mqtl$NCON[1],
            snp = merged_mqtl$SNP
          ),
          dataset2 = list(
            beta = merged_mqtl$BETA_mqtl,
            varbeta = merged_mqtl$SE_mqtl^2,
            type = "quant",
            N = 166,
            MAF = merged_mqtl$MAF,
            snp = merged_mqtl$SNP
          )
        )
        
        results_mqtl[[i]] <- data.frame(
          lead_snp = lead_snp,
          CHR = locus_chr,
          POS = lead_snps$POS[i],
          nsnps = result$summary["nsnps"],
          PP0 = result$summary["PP.H0.abf"],
          PP1 = result$summary["PP.H1.abf"],
          PP2 = result$summary["PP.H2.abf"],
          PP3 = result$summary["PP.H3.abf"],
          PP4 = result$summary["PP.H4.abf"]
        )
        
      }, error = function(e) {
        cat("mQTL error at locus", i, ":", conditionMessage(e), "\n")
      })
    }
  }
  if(nrow(gwas_locus) >= 100 & nrow(eqtl_locus) >= 100) {
    merged_eqtl <- inner_join(gwas_locus, eqtl_locus, by = c("CHR", "POS"))
    merged_eqtl <- merged_eqtl %>%
      group_by(SNP) %>%
      slice_min(P_eqtl, n = 1) %>%
      ungroup()
    
    if(nrow(merged_eqtl) >= 100) {
      
      tryCatch({
        
        result <- coloc.abf(
          dataset1 = list(
            beta = merged_eqtl$BETA_gwas,
            varbeta = merged_eqtl$SE_gwas^2,
            type = "cc",
            s = merged_eqtl$NCAS[1] / (merged_eqtl$NCAS[1] + merged_eqtl$NCON[1]),
            N = merged_eqtl$NCAS[1] + merged_eqtl$NCON[1],
            snp = merged_eqtl$SNP
          ),
          dataset2 = list(
            beta = merged_eqtl$BETA_eqtl,
            varbeta = merged_eqtl$SE_eqtl^2,
            type = "quant",
            N = 1387,
            MAF = merged_eqtl$MAF,
            snp = merged_eqtl$SNP
          )
        )
        
        results_eqtl[[i]] <- data.frame(
          lead_snp = lead_snp,
          CHR = locus_chr,
          POS = lead_snps$POS[i],
          nsnps = result$summary["nsnps"],
          PP0 = result$summary["PP.H0.abf"],
          PP1 = result$summary["PP.H1.abf"],
          PP2 = result$summary["PP.H2.abf"],
          PP3 = result$summary["PP.H3.abf"],
          PP4 = result$summary["PP.H4.abf"]
        )
        
      }, error = function(e) {
        cat("eQTL error at locus", i, ":", conditionMessage(e), "\n")
      })
    }
  }
}

#Saving the result 

mqtl_results <- bind_rows(results_mqtl)
eqtl_results <- bind_rows(results_eqtl)

#Filtering the hits 

mqtl_hits <- mqtl_results %>%
  filter(PP4 > 0.8) %>%
  arrange(desc(PP4))

eqtl_hits <- eqtl_results %>%
  filter(PP4 > 0.8) %>%
  arrange(desc(PP4))

cat("Total mQTL loci tested:", nrow(mqtl_results), "\n")
cat("Total eQTL loci tested:", nrow(eqtl_results), "\n")

write.csv(mqtl_results, "coloc_mqtl_all_results.csv", row.names = FALSE)
write.csv(eqtl_results, "coloc_eqtl_all_results.csv", row.names = FALSE)
write.csv(mqtl_hits, "coloc_mqtl_hits.csv", row.names = FALSE)
write.csv(eqtl_hits, "coloc_eqtl_hits.csv", row.names = FALSE)

#Finding hits 

print(mqtl_hits)

print(eqtl_hits)

#Mapping the genes 

mqtl_hits_genes <- mqtl_hits %>%
  left_join(
    mqtl %>% 
      select(CHR, POS, GeneAnnotation, ProbeID) %>% 
      distinct(),
    by = c("CHR", "POS")
  )
print(mqtl_hits_genes %>% select(lead_snp, CHR, POS, PP4, GeneAnnotation, ProbeID))


eqtl_hits_genes <- eqtl_hits %>%
  left_join(
    eqtl %>% 
      select(CHR, POS, gene_id) %>% 
      distinct(),
    by = c("CHR", "POS")
  )

print(eqtl_hits_genes %>% select(lead_snp, CHR, POS, PP4, gene_id))

#Null result 

#Looking for 1 mb window 

mqtl_hits_genes <- mqtl_hits %>%
  rowwise() %>%
  mutate(
    GeneAnnotation = {
      nearby <- mqtl %>%
        filter(CHR == CHR,
               POS >= POS - 500000,
               POS <= POS + 500000,
               !is.na(GeneAnnotation),
               GeneAnnotation != "") %>%
        arrange(P_mqtl) %>%
        slice(1)
      if(nrow(nearby) > 0) nearby$GeneAnnotation else NA
    }
  ) %>%
  ungroup()

print(mqtl_hits_genes %>% select(lead_snp, CHR, POS, PP4, GeneAnnotation))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = c("ENSG00000114054", "ENSG00000065609"),
  mart = mart
)
print(genes)


# Fixed version

mqtl_hits_genes <- mqtl_hits %>%
  mutate(GeneAnnotation = mapply(function(chr, pos) {
    nearby <- mqtl %>%
      filter(CHR == chr,
             POS >= pos - 500000,
             POS <= pos + 500000,
             !is.na(GeneAnnotation),
             GeneAnnotation != "") %>%
      arrange(P_mqtl) %>%
      slice(1)
    if(nrow(nearby) > 0) nearby$GeneAnnotation else NA
  }, CHR, POS))

print(mqtl_hits_genes %>% select(lead_snp, CHR, POS, PP4, GeneAnnotation))
