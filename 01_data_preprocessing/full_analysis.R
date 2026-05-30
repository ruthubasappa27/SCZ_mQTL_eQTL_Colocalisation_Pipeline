#ommmmmm
#Author: Ruthushree Basappa 
#Title: Leveraging methylation Quantitative trait loci to refine GWAS risk signals and Prioritise drug targets in Neurodevelopmental Disorders
#Running full COLOC analysis 
#Reference: Spargo et al. eLife 2023


#Installing packages 

install.packages("coloc")

install.packages("susieR")

install.packages("data.table")

install.packages("dplyr")

install.packages("ggplot2")

install.packages("tidyr")

install.packages("BEDMatrix")

#Installing Biomanagers 

if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("biomaRt")

BiocManager::install("GenomicRanges")

#Loading the downloaded packages 

library(data.table)

library(dplyr)

library(coloc)

library(susieR)

library(biomaRt)

library(ggplot2)

library(tidyr)

library(BEDMatrix)

#Checking the column names 

#FOR GWAS

gwas <- fread("~/Desktop/Thesis/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv", 
              skip = "CHROM",
              nrows = 5)
colnames(gwas)

#For mQTL

mqtl <- fread("~/Desktop/Thesis/All_Imputed_BonfSignificant_mQTLs.csv", nrows=5)

colnames(mqtl)

head(mqtl)


#Preprocessing of mQTL 

bim <- fread("~/Desktop/Thesis/COLOC-reporter/ld_reference/EUR_phase3_chr1.bim",
             col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))

head(bim)

#Reading bim file for all 22 chromosones 

bim_all <- rbindlist(lapply(1:22, function(chr) {
  fread(paste0("~/Desktop/Thesis/COLOC-reporter/ld_reference/EUR_phase3_chr", chr, ".bim"),
        col.names = c("CHR", "SNP", "CM", "POS", "A1_ref", "A2_ref"))
}))

cat("Total SNPs in reference panel:", nrow(bim_all), "\n")

#TOTAL SNPs in reference panel: 22132657

#adding rsID to the mQTL dataset 

mqtl_full <- fread("~/Desktop/Thesis/All_Imputed_BonfSignificant_mQTLs.csv")

#Remaning

mqtl_full <- mqtl_full %>% rename(
  CHR = SNP_Chr,
  POS = SNP_BP,
  EA_mqtl = SNP_Allele,
  BETA_mqtl = beta,
  P_mqtl = p.value
)

#Merging with bim file to get rsID

mqtl_full <- inner_join(mqtl_full, 
                        bim_all %>% select(CHR, POS, SNP, A1_ref, A2_ref),
                        by = c("CHR", "POS"))

#Detaching biomaRt to fix the error 

detach("package:biomaRt", unload = TRUE)

library(dplyr)

#Rerunning after the fix

mqtl_full <- inner_join(mqtl_full, 
                        bim_all %>% dplyr::select(CHR, POS, SNP, A1_ref, A2_ref),
                        by = c("CHR", "POS"))
cat("mQTL SNPs after adding rsIDs:", nrow(mqtl_full), "\n")

#Result: mQTL SNPs after adding rsIDs: 572809

cat("mQTL SNPs lost (no rsID match):", nrow(fread("~/Desktop/Thesis/All_Imputed_BonfSignificant_mQTLs.csv")) - nrow(mqtl_full), "\n")

#Result: mQTL SNPs lost (no rsID match): 85683 

# Calculate SE from beta and p-value

mqtl_full <- mqtl_full %>%
  mutate(
    Z_mqtl = qnorm(P_mqtl / 2, lower.tail = FALSE),
    SE_mqtl = abs(BETA_mqtl) / Z_mqtl
  )
# Remove any infinite or zero SE values

mqtl_full <- mqtl_full %>%
  filter(is.finite(SE_mqtl),
         SE_mqtl > 0)

# Remove duplicated SNPs 

mqtl_full <- mqtl_full %>%
  group_by(SNP) %>%
  slice_min(P_mqtl, n = 1) %>%
  ungroup()
cat("Final mQTL SNPs after QC:", nrow(mqtl_full), "\n")

#Result: Final mQTL SNP after QC: 312412

#Save the result 1

fwrite(mqtl_full, "~/Desktop/Thesis/mqtl_preprocessed.csv")

cat("Saved to mqtl_preprocessed.csv\n")


#QC for GWAS PGC3

#Reference: Winkler et al. 2014

#Reloading required library

library(data.table)

library(dplyr)

library(ggplot2)


#Loading the gwas dataset 

gwas <- fread("~/Desktop/Thesis/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",
              skip = "CHROM")

cat("Raw GWAS SNPs loaded:", nrow(gwas), "\n")

cat("Columns:", colnames(gwas), "\n")

#Result: Raw GWAS SNPs loaded: 7659767 

#Renaming to match the format 

gwas <- gwas %>% rename(
  CHR = CHROM,
  SNP = ID,
  EA = A1,
  OA = A2,
  BETA_gwas = BETA,
  SE_gwas = SE,
  P_gwas = PVAL
)

cat("Columns renamed successfully\n")

colnames(gwas)

#Checking for missingness 


cat("Missing SNP IDs:", sum(is.na(gwas$SNP)), "\n")

cat("Missing BETA:", sum(is.na(gwas$BETA_gwas)), "\n")

cat("Missing SE:", sum(is.na(gwas$SE_gwas)), "\n")

cat("Missing PVAL:", sum(is.na(gwas$P_gwas)), "\n")

cat("Missing EA:", sum(is.na(gwas$EA)), "\n")

cat("Missing OA:", sum(is.na(gwas$OA)), "\n")

cat("Missing IMPINFO:", sum(is.na(gwas$IMPINFO)), "\n")

#Imputation quality filter

gwas <- gwas %>% filter(IMPINFO >= 0.8)

cat("SNPs after IMPINFO filter:", nrow(gwas), "\n")

#Result: SNPs after IMPINFO filter: 7121961 

#MAF 

gwas <- gwas %>%
  mutate(MAF = pmin(FCAS, 1 - FCAS)) %>%
  filter(MAF >= 0.01)
cat("SNPs after MAF filter:", nrow(gwas), "\n")

#Result: SNPs after MAF filter: 7121961

#Remove Ambiguous SNPs (Palindromic)

gwas <- gwas %>%
  filter(!(EA == "A" & OA == "T"),
         !(EA == "T" & OA == "A"),
         !(EA == "C" & OA == "G"),
         !(EA == "G" & OA == "C"))

cat("SNPs after removing palindromic:", nrow(gwas), "\n")

#Result: SNPs after removing palindromic: 6041433

#Removing Duplicated values 

gwas <- gwas %>%
  filter(!duplicated(SNP))

cat("SNPs after removing duplicates:", nrow(gwas), "\n")

#Results: SNPs after removing duplicates: 6041433 

#Removing extreme values 

gwas <- gwas %>%
  filter(SE_gwas > 0,
         abs(BETA_gwas) < 10)

cat("SNPs after removing extreme values:", nrow(gwas), "\n")

#Results: SNPs after removing extreme values: 6041433 

#Saving GWAS Results 

fwrite(gwas, "~/Desktop/Thesis/gwas_qc.csv")

cat("Clean GWAS saved:", nrow(gwas), "SNPs\n")

#Filtering to genome wide significant SNPs 

gwas_sig <- gwas %>% filter(P_gwas < 5e-8)

cat("Genome-wide significant SNPs:", nrow(gwas_sig), "\n")

#Result: Genome-wide significant SNPs: 17201 

#Definining Independent Loci 

gwas_sig <- gwas_sig %>% arrange(P_gwas)
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

cat("Independent loci defined:", nrow(lead_snps), "\n")


# Create regions file in reference to Spargo pipeline 

regions <- lead_snps %>%
  mutate(
    start = POS - 500000,
    end = POS + 500000,
    start = ifelse(start < 0, 0, start)
  ) %>%
  dplyr::select(SNP, CHR, start, end)

head(regions)

cat("Total regions to test:", nrow(regions), "\n")


# Create set.regions.txt in the Spargo format 
# Format: traits(tab)region
# Region format: chromosome,start,end

regions_spargo <- regions %>%
  mutate(
    traits = "SCZ,mQTL",
    region = paste(CHR, start, end, sep = ",")
  ) %>%
  dplyr::select(traits, region)

# Adding header

write.table(
  regions_spargo,
  "~/Desktop/Thesis/COLOC-reporter/scripts/set.regions.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("set.regions.txt created with", nrow(regions_spargo), "regions\n")

# Checking first few lines

head(regions_spargo)

#Adding sample size to mQTL data 

mqtl_full <- mqtl_full %>%
  mutate(N_mqtl = 166)

# Save updated mQTL

fwrite(mqtl_full, "~/Desktop/Thesis/mqtl_preprocessed.csv")

cat("mQTL updated with N column\n")

colnames(mqtl_full)

#Since MQTL doesnt have MAF 
# Getting MAF from GWAS for overlapping SNPs

gwas_maf <- gwas %>% 
  dplyr::select(SNP, MAF)

# Adding MAF to mQTL by matching SNP ID

mqtl_full <- mqtl_full %>%
  left_join(gwas_maf, by = "SNP")

# Checking for MAF

cat("mQTL SNPs with MAF:", sum(!is.na(mqtl_full$MAF)), "\n")

cat("mQTL SNPs without MAF:", sum(is.na(mqtl_full$MAF)), "\n")

#Results: mQTL SNPs with MAF: 220889 
#Results: mQTL SNPs without MAF: 91523 

#Saving Updated mQTL file 

fwrite(mqtl_full, "~/Desktop/Thesis/mqtl_preprocessed.csv")

cat("mQTL saved with MAF column\n")

colnames(mqtl_full)

# Creating GWAS_samples.txt for Spargo pipeline

gwas_samples <- data.frame(
  ID = c("SCZ", "mQTL"),
  type = c("cc", "quant"),
  prop = c(0.5, "NA"),
  traitSD = c("NA", "NA"),
  p_col = c("P_gwas", "P_mqtl"),
  stat_col = c("BETA_gwas", "BETA_mqtl"),
  N_col = c("NEFF", "N_mqtl"),
  chr_col = c("CHR", "CHR"),
  pos_col = c("POS", "POS"),
  se_col = c("SE_gwas", "SE_mqtl"),
  snp_col = c("SNP", "SNP"),
  A1_col = c("EA", "EA_mqtl"),
  A2_col = c("OA", "A2_ref"),
  freq_col = c("MAF", "MAF"),
  traitLabel = c("Schizophrenia", "Fetal_Brain_mQTL"),
  FILEPATH = c(
    "~/Desktop/Thesis/gwas_qc.csv",
    "~/Desktop/Thesis/mqtl_preprocessed.csv"
  )
)

# Saving

write.table(
  gwas_samples,
  "~/Desktop/Thesis/COLOC-reporter/scripts/GWAS_samples.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("GWAS_samples.txt created!\n")

print(gwas_samples)

#Fixing to Run in Mac 

install.packages("optparse")

install.packages("tidyverse")

library(tidyverse)

install.packages("R.utils")

install.packages("ggrepel")

install.packages("patchwork")

install.packages("kableExtra")


#Fixing BiomarT issue 

options(BioC_mirror = "https://bioconductor.org")


#LD matrix error fix 

library(data.table)

# Load bim file for chr1

bim_chr1 <- fread("~/Desktop/Thesis/COLOC-reporter/ld_reference/EUR_phase3_chr1.bim",
                  col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))

# Load snplist used by PLINK for locus 2

snplist <- fread("~/Desktop/Thesis/COLOC-reporter/coloc/results/test_locus2_coloc/data/LDmatrix/ld_matrix.snplist",
                 header = FALSE)

cat("SNPs in PLINK snplist:", nrow(snplist), "\n")

# Check if these SNPs exist in GWAS
gwas_snps <- gwas %>% filter(CHR == 1, POS >= 28532580, POS <= 29532580)
cat("GWAS SNPs in locus:", nrow(gwas_snps), "\n")

# How many overlap?
overlap <- sum(snplist$V1 %in% gwas_snps$SNP)
cat("Overlapping SNPs:", overlap, "\n")
cat("Non-overlapping:", nrow(snplist) - overlap, "\n")

#Reading result for chromosome 1 

library(data.table)

results <- fread("~/Desktop/Thesis/COLOC-reporter/coloc/results/test_locus2_v19_coloc/tables/results_summary_coloc_abf.csv")

print(results)

#Reading result for all 201 loci 

library(data.table)

all_results <- rbindlist(lapply(
  list.files("~/Desktop/Thesis/COLOC-reporter/coloc/results", 
             pattern="results_summary_coloc_abf.csv", 
             recursive=TRUE, full.names=TRUE),
  fread), fill=TRUE)

cat("Total loci with results:", nrow(all_results), "\n")
cat("Colocalised hits PP4 > 0.8:", sum(all_results$PP.H4.abf > 0.8, na.rm=TRUE), "\n")

hits <- all_results[PP.H4.abf > 0.8][order(-PP.H4.abf)]
print(hits[, .(region, nsnps, PP.H4.abf)])

#Cleaning the duplicates

library(data.table)

library(dplyr)

all_results <- rbindlist(lapply(
  list.files("~/Desktop/Thesis/COLOC-reporter/coloc/results", 
             pattern="results_summary_coloc_abf.csv", 
             recursive=TRUE, full.names=TRUE),
  fread), fill=TRUE)

# Remove duplicates

all_results_unique <- all_results[!duplicated(region)]

cat("Total unique loci tested:", nrow(all_results_unique), "\n")

#Result: Total Unique loci tested:128

cat("Colocalised hits PP4 > 0.8:", sum(all_results_unique$PP.H4.abf > 0.8, na.rm=TRUE), "\n")

#Colocalised hits > 0.8: 6

hits <- all_results_unique[PP.H4.abf > 0.8][order(-PP.H4.abf)]

print(hits[, .(region, nsnps, PP.H4.abf)])

#Result 

#region nsnps PP.H4.abf
#<char> <int>     <num>
#  1:   Chr1:28532580-29532580   285 0.9928229
#2: Chr2:212685645-213685645    52 0.9895870
#3:   Chr6:27859632-28859632   611 0.9451987
#4:    Chr8:9532894-10532894   147 0.9131988
#5: Chr5:139833952-140833952   374 0.8811764
#6: Chr2:172456449-173456449   147 0.8160114

# Load mQTL data to get gene annotations

mqtl <- fread("~/Desktop/Thesis/mqtl_preprocessed.csv")

# For each hit find the most significant mQTL in the window

hits_genes <- hits %>%
  rowwise() %>%
  mutate(
    chr_num = as.integer(gsub("Chr(\\d+):.*", "\\1", region)),
    start = as.integer(gsub("Chr\\d+:(\\d+)-.*", "\\1", region)),
    end = as.integer(gsub("Chr\\d+:\\d+-(\\d+)", "\\1", region)),
    GeneAnnotation = {
      nearby <- mqtl %>%
        dplyr::filter(CHR == chr_num,
                      POS >= start,
                      POS <= end,
                      !is.na(GeneAnnotation),
                      GeneAnnotation != "") %>%
        arrange(P_mqtl) %>%
        slice(1)
      if(nrow(nearby) > 0) nearby$GeneAnnotation else NA
    }
  ) %>%
  ungroup()

print(hits_genes %>% dplyr::select(region, nsnps, PP.H4.abf, GeneAnnotation))


#Gene annotation attempt on exact location 

for(i in 1:nrow(hits)) {
  
  chr_num <- as.integer(gsub("Chr(\\d+):.*", "\\1", hits$region[i]))
  start <- as.integer(gsub("Chr\\d+:(\\d+)-.*", "\\1", hits$region[i]))
  end <- as.integer(gsub("Chr\\d+:\\d+-(\\d+)", "\\1", hits$region[i]))
  
  pattern <- paste0(chr_num, "_", start, "_", end, "_coloc")
  all_folders <- list.dirs("~/Desktop/Thesis/COLOC-reporter/coloc/results", 
                           recursive=FALSE)
  matched_folder <- all_folders[grepl(pattern, all_folders)]
  
  if(length(matched_folder) == 0) next
  
  snpwise_file <- file.path(matched_folder[1], 
                            "tables/coloc.abf_snpwise_PP_H4_abf.csv")
  
  if(!file.exists(snpwise_file)) next
  
  snpwise <- fread(snpwise_file)
  snpwise <- snpwise[order(-SNP.PP.H4)]
  snpwise[, cumPP := cumsum(SNP.PP.H4) / sum(SNP.PP.H4)]
  
  credible_set <- snpwise[cumPP <= 0.95]
  
  cat("=== Locus:", hits$region[i], "| PP4:", round(hits$PP.H4.abf[i], 3), "===\n")
  cat("Total SNPs:", nrow(snpwise), "\n")
  cat("95% credible set size:", nrow(credible_set), "SNPs\n")
  cat("Top SNP:", snpwise$snp[1], "| SNP.PP.H4:", round(snpwise$SNP.PP.H4[1], 4), "\n")
  
  if(snpwise$SNP.PP.H4[1] > 0.5) {
    cat("Signal: DOMINANT SPIKE - high confidence causal variant\n")
  } else if(snpwise$SNP.PP.H4[1] > 0.2) {
    cat("Signal: MODERATE - credible set approach needed\n")
  } else {
    cat("Signal: FLAT - entire credible set equally plausible\n")
  }
  cat("\n")
}


# For each locus map credible set SNPs to mQTL gene annotations
for(i in 1:nrow(hits)) {
  
  chr_num <- as.integer(gsub("Chr(\\d+):.*", "\\1", hits$region[i]))
  start <- as.integer(gsub("Chr\\d+:(\\d+)-.*", "\\1", hits$region[i]))
  end <- as.integer(gsub("Chr\\d+:\\d+-(\\d+)", "\\1", hits$region[i]))
  
  pattern <- paste0(chr_num, "_", start, "_", end, "_coloc")
  all_folders <- list.dirs("~/Desktop/Thesis/COLOC-reporter/coloc/results", 
                           recursive=FALSE)
  matched_folder <- all_folders[grepl(pattern, all_folders)]
  if(length(matched_folder) == 0) next
  
  snpwise_file <- file.path(matched_folder[1], 
                            "tables/coloc.abf_snpwise_PP_H4_abf.csv")
  if(!file.exists(snpwise_file)) next
  
  snpwise <- fread(snpwise_file)
  snpwise <- snpwise[order(-SNP.PP.H4)]
  snpwise[, cumPP := cumsum(SNP.PP.H4) / sum(SNP.PP.H4)]
  credible_set <- snpwise[cumPP <= 0.95]

  # Map credible set SNPs to mQTL genes
  
  genes <- mqtl %>%
    dplyr::filter(CHR == chr_num,
                  POS %in% credible_set$pos,
                  !is.na(GeneAnnotation),
                  GeneAnnotation != "") %>%
    dplyr::select(SNP, POS, GeneAnnotation, ProbeID, P_mqtl) %>%
    arrange(P_mqtl) %>%
    distinct(GeneAnnotation, .keep_all = TRUE)
  
  cat("=== Locus:", hits$region[i], "===\n")
  cat("Credible set SNPs:", nrow(credible_set), "\n")
  if(nrow(genes) > 0) {
    cat("Candidate genes:\n")
    print(genes %>% dplyr::select(SNP, GeneAnnotation, ProbeID))
  } else {
    cat("No exact mQTL match — will need window-based mapping\n")
  }
  cat("\n")
}

#Window mapping for CHr 6 and 8 

no_match_loci <- hits[3:4, ]

for(i in 1:nrow(no_match_loci)) {
  chr_num <- as.integer(gsub("Chr(\\d+):.*", "\\1", no_match_loci$region[i]))
  start <- as.integer(gsub("Chr\\d+:(\\d+)-.*", "\\1", no_match_loci$region[i]))
  end <- as.integer(gsub("Chr\\d+:\\d+-(\\d+)", "\\1", no_match_loci$region[i]))
  nearby_genes <- mqtl %>%
    dplyr::filter(CHR == chr_num, POS >= start, POS <= end,
                  !is.na(GeneAnnotation), GeneAnnotation != "") %>%
    arrange(P_mqtl) %>%
    distinct(GeneAnnotation, .keep_all = TRUE) %>%
    head(5)
  cat("=== Locus:", no_match_loci$region[i], "===\n")
  if(nrow(nearby_genes) > 0) {
    print(nearby_genes %>% dplyr::select(SNP, POS, GeneAnnotation, ProbeID, P_mqtl))
  } else {
    cat("No mQTL genes found in window\n")
  }
  cat("\n")
}

# For each hit load the harmonised summary stats
#  check if beta_GWAS and beta_mQTL are concordant

harmonised <- fread("~/Desktop/Thesis/COLOC-reporter/coloc/results/1_28532580_29532580_coloc/data/datasets/harmonised_sumstats.csv")

colnames(harmonised)

head(harmonised)

# Get top SNP from OPRD1 credible set

top_snp <- "rs61783570"

directionality <- harmonised %>%
  dplyr::filter(snp == top_snp) %>%
  dplyr::select(snp, trait, beta, SE, pvalues, cs)

print(directionality)

# Define top SNPs for each hit

library(data.table)
library(dplyr)

top_snps <- c("rs61783570",   
              "rs7607363",    
              "rs3118359",    
              "rs4633059",    
              "rs2530240",    
              "rs4667693")    

locus_folders <- c(
  "1_28532580_29532580",
  "2_212685645_213685645", 
  "6_27859632_28859632",
  "8_9532894_10532894",
  "5_139833952_140833952",
  "2_172456449_173456449"
)

for(i in 1:length(top_snps)) {
  harm_file <- paste0("~/Desktop/Thesis/COLOC-reporter/coloc/results/",
                      locus_folders[i], "_coloc/data/datasets/harmonised_sumstats.csv")
  
  if(!file.exists(harm_file)) {
    cat("No harmonised file for:", locus_folders[i], "\n\n")
    next
  }
  
  harm <- fread(harm_file)
  
  dir_result <- harm %>%
    dplyr::filter(snp == top_snps[i]) %>%
    dplyr::select(snp, trait, beta, pvalues, cs)
  
  scz_beta <- dir_result$beta[dir_result$trait == "SCZ"]
  mqtl_beta <- dir_result$beta[dir_result$trait == "mQTL"]
  
  cat("=== Gene", i, "| Locus:", locus_folders[i], "| SNP:", top_snps[i], "===\n")
  print(dir_result)
  
  if(length(scz_beta) > 0 && length(mqtl_beta) > 0) {
    if(sign(scz_beta) == sign(mqtl_beta)) {
      cat("Direction: CONCORDANT\n")
    } else {
      cat("Direction: DISCORDANT\n")
    }
  }
  cat("\n")
}
}

#RESULTS 

#Running SMR 

library(data.table)

library(dplyr)

# Load QC GWAS data 

gwas <- fread("~/Desktop/Thesis/gwas_qc.csv")

# Create SMR format GWAS file

gwas_smr <- gwas %>%
  dplyr::select(
    SNP = SNP,
    A1 = EA,
    A2 = OA,
    freq = MAF,
    b = BETA_gwas,
    se = SE_gwas,
    p = P_gwas,
    N = NEFF
  )

fwrite(gwas_smr, "~/Desktop/Thesis/gwas_smr.txt", sep="\t")

cat("GWAS SMR file created:", nrow(gwas_smr), "SNPs\n")

#converting mqtl into ESD and BESD format 


# Format: Chr ProbeID GeneticDistance ProbeBP Gene Strand

mqtl <- fread("~/Desktop/Thesis/mqtl_preprocessed.csv")

epi <- mqtl %>%
  dplyr::select(
    Chr = DNAm_CHR,
    ProbeID = ProbeID,
    ProbeBP = DNAm_BP,
    Gene = GeneAnnotation
  ) %>%
  distinct(ProbeID, .keep_all = TRUE) %>%
  mutate(
    GeneticDistance = 0,
    Strand = "+"
  ) %>%
  dplyr::select(Chr, ProbeID, GeneticDistance, ProbeBP, Gene, Strand)

fwrite(epi, "~/Desktop/Thesis/mqtl_smr.epi", 
       sep="\t", col.names=FALSE)
cat("EPI file created:", nrow(epi), "probes\n")

# Creating ESI file: one row per SNP
# Format: Chr SNP GeneticDistance BP A1 A2 Freq

esi <- mqtl %>%
  dplyr::select(
    Chr = CHR,
    SNP = SNP,
    BP = POS,
    A1 = EA_mqtl,
    A2 = A2_ref,
    Freq = MAF
  ) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  mutate(
    GeneticDistance = 0,
    Freq = ifelse(is.na(Freq), 0.5, Freq)
  ) %>%
  dplyr::select(Chr, SNP, GeneticDistance, BP, A1, A2, Freq)

fwrite(esi, "~/Desktop/Thesis/mqtl_smr.esi",
       sep="\t", col.names=FALSE)
cat("ESI file created:", nrow(esi), "SNPs\n")

# creating ESD file: effect sizes linking SNPs to probes
# Format: one file per probe containing SNP effects


esd <- mqtl %>%
  dplyr::select(
    ProbeID = ProbeID,
    SNP = SNP,
    A1 = EA_mqtl,
    A2 = A2_ref,
    Freq = MAF,
    b = BETA_mqtl,
    se = SE_mqtl,
    p = P_mqtl
  ) %>%
  mutate(
    Freq = ifelse(is.na(Freq), 0.5, Freq)
  )

fwrite(esd, "~/Desktop/Thesis/mqtl_smr_esd.txt",
       sep="\t", col.names=TRUE)
cat("ESD file created:", nrow(esd), "SNP-probe pairs\n")


#Using PLINK for BESD 

#Debugging 

# Creating proper SMR format file

smr_format <- mqtl %>%
  dplyr::select(
    ProbeID = ProbeID,
    Chr = DNAm_CHR,
    ProbeBP = DNAm_BP,
    Gene = GeneAnnotation,
    SNP = SNP,
    SNPChr = CHR,
    SNPPos = POS,
    A1 = EA_mqtl,
    A2 = A2_ref,
    Freq = MAF,
    b = BETA_mqtl,
    se = SE_mqtl,
    p = P_mqtl
  ) %>%
  mutate(
    Freq = ifelse(is.na(Freq), 0.5, Freq),
    Orientation = "+"
  ) %>%
  dplyr::select(ProbeID, Chr, ProbeBP, Gene, Orientation, 
                SNP, SNPChr, SNPPos, A1, A2, Freq, b, se, p)

fwrite(smr_format, "~/Desktop/Thesis/mqtl_smr_format.txt", sep="\t")
cat("SMR format file created:", nrow(smr_format), "rows\n")

#Retrying 

library(data.table)

library(dplyr)

mqtl <- fread("~/Desktop/Thesis/mqtl_preprocessed.csv")

esi_chr1 <- mqtl %>%
  dplyr::filter(CHR == 1) %>%
  dplyr::select(CHR, SNP, POS, EA_mqtl, A2_ref, MAF) %>%
  distinct(SNP, .keep_all=TRUE) %>%
  mutate(
    GeneticDistance = 0,
    MAF = ifelse(is.na(MAF), 0.5, MAF)
  ) %>%
  dplyr::select(CHR, SNP, GeneticDistance, POS, EA_mqtl, A2_ref, MAF)

fwrite(esi_chr1, "~/Desktop/Thesis/test_chr1.esi", 
       sep="\t", col.names=FALSE)
cat("Chr1 ESI:", nrow(esi_chr1), "SNPs\n")

library(data.table)

library(dplyr)

# Read finemapping results for all 6 hit loci

locus_folders <- c(
  "1_28532580_29532580",
  "2_212685645_213685645", 
  "6_27859632_28859632",
  "8_9532894_10532894",
  "5_139833952_140833952",
  "2_172456449_173456449"
)

for(locus in locus_folders) {
  finemapping_file <- paste0("~/Desktop/Thesis/COLOC-reporter/coloc/results/",
                             locus, "_coloc/tables/results_summary_finemapping.csv")
  
  if(file.exists(finemapping_file)) {
    fm <- fread(finemapping_file)
    cat("=== Locus:", locus, "===\n")
    print(fm %>% dplyr::select(trait, cs, pip_maxSNP, 
                               pip_maxSNP_nearestDownstreamGene,
                               pip_maxSNP_nearestUpstreamGene))
    cat("\n")
  }
}

library(biomaRt)

library(dplyr)


ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Top SNPs for each hit
top_snps <- data.frame(
  gene = c("OPRD1", "ERBB4", "ZNF389", "Chr8", "WDR55", "DYNC1I2"),
  snp = c("rs61783570", "rs7607363", "rs3118359", 
          "rs4633059", "rs2530240", "rs4667693"),
  chr = c(1, 2, 6, 8, 5, 2),
  pos = c(29150307, 213402705, 28751727, 
          9932468, 140024042, 172530949)
)

# Get gene information for each SNP position
for(i in 1:nrow(top_snps)) {
  cat("=== Gene:", top_snps$gene[i], "| SNP:", top_snps$snp[i], "===\n")
  
  results <- getBM(
    attributes = c('external_gene_name', 'chromosome_name', 
                   'start_position', 'end_position', 'gene_biotype'),
    filters = c('chromosome_name', 'start', 'end'),
    values = list(top_snps$chr[i], 
                  top_snps$pos[i] - 10000, 
                  top_snps$pos[i] + 10000),
    mart = ensembl
  )
  
  print(results)
  cat("\n")
  
  
  library(httr)
  
  library(jsonlite)
  
  library(dplyr)
  
  # Function to get SNP annotation from myvariant.info
  
  get_snp_annotation <- function(rsid) {
    url <- paste0("https://myvariant.info/v1/variant/", rsid, 
                  "?assembly=hg19&fields=cadd.gene,cadd.consequence,dbsnp.gene")
    
    response <- GET(url)
    
    if(status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding="UTF-8"))
      return(data)
    } else {
      return(NULL)
    }
  }
  
  # Top SNPs
  
  snps <- c("rs61783570", "rs7607363", "rs3118359", 
            "rs4633059", "rs2530240", "rs4667693")
  genes <- c("OPRD1", "ERBB4", "ZNF389", "Chr8", "WDR55", "DYNC1I2")
  
  for(i in 1:length(snps)) {
    cat("=== Candidate gene:", genes[i], "| SNP:", snps[i], "===\n")
    result <- get_snp_annotation(snps[i])
    
    if(!is.null(result)) {
      
      # Extract gene and consequence
      
      tryCatch({
        if(!is.null(result$cadd$gene$genename)) {
          cat("Gene:", result$cadd$gene$genename, "\n")
        }
        if(!is.null(result$cadd$consequence)) {
          cat("Consequence:", result$cadd$consequence, "\n")
        }
        if(!is.null(result$dbsnp$gene$name)) {
          cat("dbSNP gene:", result$dbsnp$gene$name, "\n")
        }
      }, error = function(e) {
        cat("Could not parse annotation\n")
      })
    }
    cat("\n")
  }
}

library(data.table)
library(dplyr)

mqtl <- fread("~/Desktop/Thesis/mqtl_preprocessed.csv")

# Check what the mQTL data says about rs2530240

snp_info <- mqtl %>%
  dplyr::filter(SNP == "rs2530240") %>%
  dplyr::select(SNP, CHR, POS, ProbeID, DNAm_CHR, DNAm_BP, GeneAnnotation)

print(snp_info)

# Also check all SNPs in the WDR55 credible set

# Load the snpwise file for Chr5 locus

snpwise <- fread("~/Desktop/Thesis/COLOC-reporter/coloc/results/5_139833952_140833952_coloc/tables/coloc.abf_snpwise_PP_H4_abf.csv")

# Get credible set SNPs

snpwise <- snpwise[order(-SNP.PP.H4)]
snpwise[, cumPP := cumsum(SNP.PP.H4) / sum(SNP.PP.H4)]
credible_set <- snpwise[cumPP <= 0.95]

cat("Credible set SNPs for Chr5 locus:\n")

print(credible_set %>% dplyr::select(snp, pos, SNP.PP.H4))


# Check the exact position of the CpG relative to WDR55

cat("CpG probe cg26395211 position: chr5:", 140044315, "\n")

cat("SNP rs2530240 position: chr5:", 140060547, "\n")

cat("Distance SNP to CpG:", abs(140060547 - 140044315), "bp\n")

# Now check via myvariant for the CpG position

library(httr)

library(jsonlite)

# Query the CpG position

url <- "https://myvariant.info/v1/query?q=chr5:140044315&assembly=hg19&fields=cadd.gene"
response <- GET(url)
data <- fromJSON(content(response, "text", encoding="UTF-8"))
print(data$hits)

library(httr)
library(jsonlite)

# Query rs2530240 explicitly on hg19

url <- "https://myvariant.info/v1/variant/rs2530240?assembly=hg19&fields=cadd,dbsnp,hg19"

response <- GET(url)
data <- fromJSON(content(response, "text", encoding="UTF-8"))

# Check the actual coordinates returned

cat("hg19 position:", data$hg19$start, "\n")

cat("Gene:", data$cadd$gene$genename, "\n")

cat("Consequence:", data$cadd$consequence, "\n")
