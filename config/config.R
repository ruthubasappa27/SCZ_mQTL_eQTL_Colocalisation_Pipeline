# config.R
# Paths and analysis settings used by all scripts.

suppressMessages(library(here))
here::i_am("config/config.R")

PATHS <- list(
  # Raw downloads
  gwas_raw  = here("data", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv"),
  mqtl_raw  = here("data", "All_Imputed_BonfSignificant_mQTLs.csv"),
  eqtl_raw  = here("data", "all_eqtls_gene.txt"),
  eqtl_snp_pos = here("data", "snp_positions.bed"),
  # CoDEx exports (codex_mqtl.csv, codex_eqtl.csv)
  codex_dir = here("data", "codex"),

  # 1000 Genomes Phase 3 EUR, per chromosome: EUR_phase3_chr1.bed/.bim/.fam ... chr22
  ld_prefix = here("references", "ld_reference", "EUR_phase3_chr"),
  # Hannon fetal brain mQTL in SMR format (FB_Brain_2.besd/.epi/.esi)
  besd_mqtl = here("references", "FB_Brain_2", "FB_Brain_2"),
  smr_ref   = here("references", "smr_ld_chrbp"),
  # downloaded on first use: UCSC hg19 refGene; Kang et al. (2011) GSE25219 and GPL5175
  annotation_dir = here("references", "annotation"),
  kang_dir  = here("references", "kang2011"),

  # Programs
  plink = here("tools", "plink", "plink"),
  smr   = here("tools", "smr", "smr"),

  # Outputs
  processed = here("data", "processed"),
  results   = here("results"),
  logs      = here("logs"),
  tmp       = here("tmp")
)

SETTINGS <- list(
  # GWAS QC
  gwas_info_min    = 0.8,     # imputation INFO
  maf_min          = 0.01,    # MAF from cases+controls combined
  drop_palindromic = TRUE,    # remove A/T and C/G SNPs (strand cannot be checked)
  max_abs_beta     = 10,
  gwas_n_cases     = 53386,
  gwas_n_controls  = 77258,
  # NEFF = total effective sample size (PGC3 file header); case fraction 0.5 for effective N
  gwas_neff_multiplier = 1,
  gwas_case_fraction   = 0.5,

  # Locus definition
  gwas_p_sig        = 5e-8,
  locus_half_window = 500000,  # lead SNP +/- 500 kb
  mhc_chr   = 6,
  mhc_start = 25000000,        # extended MHC, hg19
  mhc_end   = 34000000,
  mhc_action = "exclude",      # "exclude" or "keep" (MHC analysed as one locus)

  # QTL
  n_mqtl = 166,
  n_eqtl = 120,
  qtl_feature_min_p = 1e-4,    # feature tested if its smallest QTL p in the locus is below this
  min_overlap_snps  = 10,      # minimum SNPs shared by GWAS, QTL and LD panel
  max_freq_diff     = 0.2,     # SNP removed if its allele frequency differs from the LD panel by more (as in SMR)

  # Colocalisation (same priors for coloc.abf and coloc.susie)
  p1  = 1e-4,
  p2  = 1e-4,
  p12 = 1e-5,
  pp4_threshold      = 0.8,
  susie_max_signals  = 10,
  # credible-set purity computed over all SNPs in the set (susieR default subsamples 100)
  susie_n_purity     = Inf,
  seed               = 20261023,  # set before every SuSiE fit
  credible_coverage  = 0.95,
  sensitivity_p12    = c(1e-6, 5e-6, 1e-5, 5e-5),

  # SMR / HEIDI
  smr_peqtl_main        = 5e-8,  # standard SMR instrument threshold (primary)
  smr_peqtl_sensitivity = 1e-4,  # sensitivity analysis
  smr_threads           = 4
)
