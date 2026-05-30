# SCZ mQTL Drug Target Prioritisation Pipeline

## Overview
This repository contains the computational pipeline for identifying and prioritising drug targets in schizophrenia (SCZ) through integration of genome-wide association study (GWAS) summary statistics with fetal brain methylation quantitative trait loci (mQTL) data.

**Author:** Ruthushree P Basappa  
**Institution:** University of Galway  
**Supervisor:** Prof. Derek Morris  
**Degree:** MSc Clinical Neuroscience  

## Study Design
- **GWAS data:** PGC3 SCZ GWAS (Trubetskoy et al. 2022, Nature) — n=127,906
- **mQTL data:** Hannon et al. fetal brain mQTL (2016, Nature Neuroscience) — n=166
- **Genome build:** GRCh37/hg19
- **Pipeline:** COLOC-reporter (Spargo et al. 2023, eLife)

## Key Results
6 colocalised loci identified (PP4 >= 0.8):

| Gene | PP4 | p_SMR | p_HEIDI |
|------|-----|-------|---------|
| OPRD1 | 0.993 | 1.19e-05 | NA |
| ERBB4 | 0.990 | 5.61e-07 | 0.108 |
| ZNF389 | 0.945 | NA | NA |
| PIK3C2B | 0.913 | 1.13e-04 | 0.285 |
| WDR55 | 0.881 | 3.51e-05 | 0.604 |
| DYNC1I2 | 0.816 | 1.42e-04 | 0.825 |

## Repository Structure
## Requirements

### Software
- R version 4.5.1
- PLINK v1.90
- SMR v1.03

### R Packages
- data.table
- dplyr
- coloc (v5.2.3)
- susieR
- httr
- jsonlite

## Data Sources
All datasets are publicly available:

1. **PGC3 SCZ GWAS:** https://pgc.unc.edu/for-researchers/download-results/
2. **Hannon fetal brain mQTL:** http://epigenetics.essex.ac.uk/mQTL/
3. **Hannon BESD format mQTL:** https://yanglab.westlake.edu.cn/software/smr/#DataResource
4. **1000 Genomes LD reference:** https://yanglab.westlake.edu.cn/software/smr/#DataResource
5. **COLOC-reporter pipeline:** https://github.com/ThomasSpargo/COLOC-reporter

## Pipeline Steps

### Step 1 — Data Preprocessing
Run `01_data_preprocessing/full_analysis.R`
- GWAS QC: removes SNPs with IMPINFO < 0.8, MAF < 0.01, palindromic SNPs
- mQTL preprocessing: standardises columns, calculates SE, assigns rsIDs
- Output: gwas_qc.csv (6,041,433 SNPs), mqtl_preprocessed.csv (312,412 SNPs)

### Step 2 — Locus Definition
Run within `01_data_preprocessing/full_analysis.R`
- Distance-based greedy algorithm, 1Mb windows
- Output: 201 independent loci

### Step 3 — Colocalisation
Run `03_colocalisation/run_coloc_all.sh`
- COLOC-reporter pipeline with trySusie mode
- PP4 >= 0.8 threshold
- Output: results in COLOC-reporter/coloc/results/

### Step 4 — Results Analysis
Run `04_results_analysis/Project Sample run 1.R`
- Credible set construction
- Directionality analysis
- Variant annotation

### Step 5 — Functional Validation
- CADD scores: myvariant.info API (hg19)
- ClinVar: myvariant.info API
- Cell type specificity: Polioudakis et al. 2019 (UCLA Geschwind lab)
- LD analysis: LDlink (ldlink.nci.nih.gov)

### Step 6 — SMR Analysis
Run SMR using:
```bash
./smr_x86 \
  --bfile g1000_eur_chrbp/g1000_eur \
  --gwas-summary gwas_smr_chrbp.txt \
  --beqtl-summary FB_Brain_2 \
  --out smr_all_genes \
  --extract-probe all_probes.txt \
  --peqtl-smr 5e-8 \
  --diff-freq 0.4 \
  --diff-freq-prop 0.6 \
  --thread-num 4
```

## References
1. Trubetskoy V, et al. Nature. 2022;604:502-508.
2. Hannon E, et al. Nature Neuroscience. 2016;19:1568-1579.
3. Spargo TP, et al. eLife. 2023;12:e88768.
4. Giambartolomei C, et al. PLOS Genetics. 2014;10:e1004383.
5. Wallace C. PLOS Genetics. 2021;17:e1009440.
6. Zhu Z, et al. Nature Genetics. 2016;48:481-487.
7. Wang G, et al. J Royal Statistical Society B. 2020;82:1273-1300.
8. 1000 Genomes Project Consortium. Nature. 2015;526:68-74.
9. Polioudakis D, et al. Neuron. 2019;103:785-801.
10. Machiela MJ, Chanock SJ. Bioinformatics. 2015;31:3555-3557.
11. Rentzsch P, et al. Nucleic Acids Research. 2019;47:D886-D894.
