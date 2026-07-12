# SCZ mQTL eQTL Colocalisation Pipeline

## Overview

This repository contains the computational pipeline used to prioritise genes for schizophrenia (SCZ) by integrating genome-wide association study (GWAS) summary statistics with fetal brain methylation quantitative trait loci (mQTL) and expression quantitative trait loci (eQTL) data.

The project uses SCZ–fetal brain colocalisation as the primary statistical analysis, followed by causal-support analysis (SMR/HEIDI). The workflow is designed to produce a conservative, evidence-based, genetically prioritised set of genes for future functional investigation.

**Author:** Ruthushree P Basappa  
**Institution:** University of Galway  
**Supervisor:** Prof. Derek Morris  
**Degree:** MSc Clinical Neuroscience  

## Study design

- **GWAS:** PGC3 schizophrenia GWAS (Trubetskoy et al. 2022), 53,386 cases; 77,258 controls
- **mQTL:** Hannon et al. fetal brain mQTL (2016), n=166
- **eQTL:** O'Brien et al. fetal brain eQTL (2018), n=120
- **Genome build:** GRCh37 / hg19
- **Primary colocalisation framework:** COLOC-reporter (Spargo et al. 2024)

## Key results

### mQTL colocalisation

Six loci showed strong evidence of colocalisation between SCZ GWAS and fetal brain mQTL signals using a PP4 threshold of 0.8.

| Gene | PP4 | p_SMR | p_HEIDI |
|------|-----:|------:|--------:|
| OPRD1 | 0.993 | 1.19e-05 | NA |
| ERBB4 | 0.990 | 5.61e-07 | 0.108 |
| ZNF389 | 0.945 | NA | NA |
| PIK3C2B | 0.913 | 1.13e-04 | 0.285 |
| WDR55 | 0.881 | 3.51e-05 | 0.604 |
| DYNC1I2 | 0.816 | 1.42e-04 | 0.825 |

### eQTL colocalisation

Locus-level colocalisation was tested across all 201 independent GWAS loci; four loci passed the PP4 ≥ 0.8 threshold. Per-gene colocalisation was then performed at these four loci (following Dobbyn et al. 2018), identifying seven gene-level signals across three loci.

| Gene | PP4 | p_SMR | p_HEIDI |
|------|-----:|------:|--------:|
| CYP2D6 | 0.987 | 7.36e-05 | 2.03e-04 |
| NAGA | 0.965 | 2.79e-04 | 0.888 |
| SMDT1 | 0.955 | 1.02e-04 | 0.0647 |
| WBP2NL | 0.887 | 4.76e-05 | 0.0176 |
| HSPA9 | 0.881 | 4.39e-04 | 0.416 |
| ENSG00000227370 | 0.850 | NA | NA |
| KRT18P46 | 0.803 | 1.61e-03 | 0.375 |

HEIDI was used as a heterogeneity check in the SMR framework. p_HEIDI > 0.05 is consistent with a single shared causal variant (pleiotropy); p_HEIDI < 0.05 suggests linkage between distinct causal variants.

Full colocalisation results for significant loci are not stored in this repository; they are available via the Data Availability statement in the associated manuscript.

## Repository structure

```text
SCZ_mQTL_eQTL_Colocalisation_Pipeline/
├── data/
│   └── README.md
├── references/
│   └── README.md
├── results/
│   └── README.md
├── 01_data_preprocessing/
│   └── preprocess_and_define_loci.R
├── 01b_data_preprocessing_eqtl/
│   └── eqtl_R_analysis.R
├── 02_colocalisation/
│   ├── GWAS_samples.txt
│   └── run_coloc_all.sh
├── 02b_colocalisation_eqtl/
│   ├── GWAS_samples_eqtl.txt
│   └── run_coloc_all_eqtl.sh
├── 03_results_analysis/
│   └── analyse_coloc_results.R
├── 03b_results_analysis_eqtl/
│   └── analyse_coloc_results_eqtl.R
├── 04_SMR_analysis/
│   ├── all_probes.txt
│   ├── prepare_gwas_smr.R
│   └── run_smr.sh
├── 04b_SMR_analysis_eqtl/
│   ├── eqtl_flist.txt
│   ├── build_eqtl_besd.sh
│   └── run_smr_eqtl.sh
└── README.md
```

## Requirements

### Software

- R 4.5.1
- PLINK v1.90
- SMR v1.03

### R packages

- data.table
- dplyr
- coloc
- susieR
- httr
- jsonlite
- here

### Operating system

- Tested on macOS Apple Silicon
- Linux-compatible with the appropriate SMR binary

## Data sources

All datasets used in this pipeline are publicly available.

1. **PGC3 SCZ GWAS:** [PGC download page](https://pgc.unc.edu/for-researchers/download-results/)
2. **Hannon fetal brain mQTL:** [Essex mQTL resource](https://epigenetics.essex.ac.uk/mQTL/) and [PubMed article](https://pubmed.ncbi.nlm.nih.gov/26619357/)
3. **Hannon BESD-format mQTL:** [SMR data resources](https://yanglab.westlake.edu.cn/software/smr/#DataResource)
4. **O'Brien fetal brain eQTL summary statistics:** [figshare dataset](https://doi.org/10.6084/m9.figshare.6881825), [Genome Biology article](https://doi.org/10.1186/s13059-018-1567-1)
5. **1000 Genomes Phase 3 EUR LD reference panel:** [PLINK format](https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe), [SMR format](https://yanglab.westlake.edu.cn/software/smr/)
6. **COLOC-reporter pipeline:** [GitHub repository](https://github.com/ThomasPSpargo/COLOC-reporter)

## Workflow

### 1. Data preprocessing

**mQTL:** Run `01_data_preprocessing/preprocess_and_define_loci.R`.
**eQTL:** Run `01b_data_preprocessing_eqtl/eqtl_R_analysis.R`.

These steps perform:
- GWAS QC, including filtering by imputation quality, MAF, and palindromic SNP removal
- mQTL/eQTL preprocessing, including column standardisation, standard error calculation, and allele harmonisation
- Export of cleaned GWAS and mQTL/eQTL summary-statistics files

**Outputs:**
- `data/gwas_qc.csv`
- `data/mqtl_preprocessed.csv`
- `data/eqtl_preprocessed.csv`

### 2. Locus definition

Performed within `01_data_preprocessing/preprocess_and_define_loci.R`.

A distance-based greedy approach was used to define independent GWAS loci using 1 Mb windows centered on lead SNPs.

**Output:** `201` independent loci

### 3. Colocalisation

**mQTL:** Run `02_colocalisation/run_coloc_all.sh`.
**eQTL:** Run `02b_colocalisation_eqtl/run_coloc_all_eqtl.sh`.

This step uses COLOC-reporter to test for shared causal variants at each locus.

- `coloc.susie` is used when SuSiE fine-mapping converges
- `coloc.abf` is used as a fallback when SuSiE does not converge
- Loci with PP4 ≥ 0.8 are considered colocalised

For eQTL, colocalisation was first tested at the locus level across all 201 loci. Loci passing PP4 ≥ 0.8 were then carried forward for per-gene colocalisation (following Dobbyn et al. 2018), run independently for each gene with available eQTL data within the locus window, to identify which specific gene(s) drove the colocalised signal. This per-gene step is performed in `03b_results_analysis_eqtl/analyse_coloc_results_eqtl.R`.

### 4. Results analysis

**mQTL:** Run `03_results_analysis/analyse_coloc_results.R`.
**eQTL:** Run `03b_results_analysis_eqtl/analyse_coloc_results_eqtl.R`.

These steps include:
- per-gene colocalisation (eQTL only)
- candidate credible set construction
- directionality analysis
- locus-level and gene-level prioritisation

### 5. SMR analysis

SMR was used as an orthogonal analysis to evaluate whether the SCZ GWAS signal and fetal brain mQTL/eQTL signal at each colocalised locus are consistent with a shared underlying variant.

**mQTL:** Run `04_SMR_analysis/prepare_gwas_smr.R`, then `04_SMR_analysis/run_smr.sh`.

**eQTL:** No pre-built BESD file was available for the O'Brien eQTL dataset (unlike FB_Brain_2 for mQTL), so per-gene BESD files were constructed directly. Run `04b_SMR_analysis_eqtl/build_eqtl_besd.sh`, then `04b_SMR_analysis_eqtl/run_smr_eqtl.sh`. Custom BESD files are referenced via the `--eqtl-flist` flag (`04b_SMR_analysis_eqtl/eqtl_flist.txt`), using an instrument significance threshold of 1×10⁻⁴ to maximise instrument availability given the smaller eQTL sample size (n=120), and the `--disable-freq-ck` flag following manual verification that frequency discordances reflected allele coding convention differences rather than genuine mismatches. SMR/HEIDI results are treated as orthogonal supporting evidence rather than a primary inclusion criterion.

The Hannon fetal brain mQTL BESD file and the O'Brien eQTL data both use chr:bp SNP identifiers rather than rsIDs. GWAS summary statistics and the LD reference panel were converted from rsID to chr:bp format before SMR analysis to ensure consistent SNP matching.

## Notes

- This repository is under active development
- The current workflow prioritises specificity and reproducibility over maximal sensitivity
- Conservative harmonisation choices may exclude some borderline variants, but retained loci represent high-confidence genes
- This repository is intended for computational, genetically-prioritised gene identification and does not claim experimental validation or therapeutic target confirmation

## References

1. Trubetskoy V, et al. Nature. 2022;604:502-508. [doi](https://doi.org/10.1038/s41586-022-04434-5)
2. Hannon E, et al. Nature Neuroscience. 2016;19:48-54. [doi](https://doi.org/10.1038/nn.4182)
3. O'Brien HE, et al. Genome Biology. 2018;19:194. [doi](https://doi.org/10.1186/s13059-018-1567-1)
4. Spargo TP, et al. eLife. 2024;12:RP88768. [doi](https://doi.org/10.7554/eLife.88768)
5. Giambartolomei C, et al. PLoS Genetics. 2014;10:e1004383. [doi](https://doi.org/10.1371/journal.pgen.1004383)
6. Wallace C. PLoS Genetics. 2021;17:e1009440. [doi](https://doi.org/10.1371/journal.pgen.1009440)
7. Zhu Z, et al. Nature Genetics. 2016;48:481-487. [doi](https://doi.org/10.1038/ng.3538)
8. Dobbyn A, et al. American Journal of Human Genetics. 2018;102:1169-1184. [doi](https://doi.org/10.1016/j.ajhg.2018.04.011)
9. 1000 Genomes Project Consortium. Nature. 2015;526:68-74. [doi](https://doi.org/10.1038/nature15393)
