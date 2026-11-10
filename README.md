# Colocalisation of schizophrenia GWAS with fetal brain mQTL and eQTL

This repository contains the code for the analysis of schizophrenia risk loci with fetal brain
DNA methylation and gene expression. The schizophrenia GWAS (Trubetskoy et al. 2022) is
colocalised with fetal brain methylation QTL (Hannon et al. 2016) and expression QTL
(O'Brien et al. 2018). Colocalised features are then tested with SMR and HEIDI (Zhu et al. 2016).
Genes annotated to the prioritised CpG sites are described by their developmental expression
(Kang et al. 2011) and cell-type expression (Polioudakis et al. 2019).

Author: Ruthushree P Basappa. Supervisor: Prof. Derek Morris. University of Galway.

## Requirements

Versions used: R 4.5.1, data.table 1.18.6.1, coloc 5.2.3, susieR 0.14.2, PLINK 1.9, SMR 1.03.

```r
install.packages(c("data.table", "here", "R.utils", "coloc", "susieR", "jsonlite", "ggplot2",
                   "gtable", "flextable", "officer", "openxlsx", "BiocManager"))
BiocManager::install(c("GEOquery", "org.Hs.eg.db", "AnnotationDbi"))
```

PLINK 1.9 goes in `tools/plink/plink` and SMR 1.03 in `tools/smr/smr`.

## Input data

The input data are not stored in this repository. All file paths are set in `config/config.R`.

| File | Source |
|---|---|
| `data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv` | PGC3 schizophrenia GWAS, European ancestry (Trubetskoy et al. 2022) |
| `data/All_Imputed_BonfSignificant_mQTLs.csv` | Fetal brain mQTL (Hannon et al. 2016) |
| `data/all_eqtls_gene.txt`, `data/snp_positions.bed` | Fetal brain eQTL and SNP positions (O'Brien et al. 2018) |
| `references/FB_Brain_2/FB_Brain_2.besd/.epi/.esi` | Fetal brain mQTL in SMR format (Hannon et al. 2016) |
| `references/ld_reference/EUR_phase3_chr1-22.bed/.bim/.fam` | 1000 Genomes Phase 3, European samples |
| `data/codex/codex_mqtl.csv`, `codex_eqtl.csv` | CoDEx summary exports for the prioritised genes (Polioudakis et al. 2019) |

These files are downloaded by the scripts on first use: UCSC hg19 refGene (step 8),
GEO GSE25219 and GPL5175 (step 12). Steps 8 and 10 query the Ensembl GRCh37 REST API.

## How to run

Run from the project folder. `run_all.sh` runs every step in order.

| Step | Script | What it does |
|---|---|---|
| 0 | `00_check_setup.R` | Checks R packages, programs and input files |
| 1 | `01_qc_gwas.R` | GWAS QC and alignment to the 1000 Genomes reference panel |
| 2 | `02_define_loci.R` | Loci: lead SNP +/- 500 kb, overlapping windows merged, MHC excluded |
| 3 | `03_qc_mqtl.R` | mQTL QC and alignment to the reference panel |
| 4 | `04_qc_eqtl.R` | eQTL QC and alignment to the reference panel (uses the loci from step 2) |
| 5 | `05_colocalisation.R mqtl` / `eqtl` | coloc.abf and coloc.susie for every locus and feature; prior sensitivity |
| 6 | `06_results.R mqtl` / `eqtl` | Locus overview, colocalised features, credible sets, direction of effect |
| 7 | `07_smr.R mqtl` / `eqtl` | SMR and HEIDI at p < 5e-8 (primary) and p < 1e-4 (sensitivity) |
| 8 | `08_annotate_genes.R` | Gene annotation of the colocalised features |
| 9 | `09_audit_direction.R` | Checks each direction of effect against the original files |
| 10 | `10_verify_rsids.R` | Checks every reported SNP against Ensembl GRCh37 |
| 11 | `11_prioritise.R` | Prioritised features and their genes |
| 12 | `12_developmental_expression.R` | Developmental expression (Kang et al. 2011) |
| 13 | `13_cell_type_expression.R` | Cell-type expression (CoDEx) |
| 14 | `14_main_tables.R` | Table 2 |
| 15 | `15_workflow_diagram.R` | Figure 1 |
| 16 | `16_flow_diagram.R` | Supplementary Figure S1 |
| 17 | `17_main_figures.R` | Figures 2 and 3 |
| 18 | `18_supplementary_workbook.R` | Supplementary Tables S1-S9 |

All settings (QC thresholds, priors, PP4 threshold, SMR thresholds) are in `config/config.R`.
Each script writes a log to `logs/` and its counts to `results/qc_counts.tsv`.

## Outputs

| Manuscript item | File |
|---|---|
| Table 2 | `results/tables/Table2.docx` |
| Figure 1 | `results/figures/Figure1.pdf` |
| Figures 2 and 3 | `results/figures/Figure2.pdf`, `Figure3.pdf` |
| Supplementary Figure S1 | `results/figures/Supplementary_Figure_S1.pdf` |
| Supplementary Tables S1-S9 | `results/supplementary/Supplementary_Tables.xlsx` |

Intermediate results are in `results/loci/`, `results/coloc_mqtl/`, `results/coloc_eqtl/`,
`results/smr_mqtl/`, `results/smr_eqtl/`, `results/rsid_check/` and `results/characterisation/`.

## References

- Giambartolomei C, et al. (2014). *PLoS Genetics*, 10, e1004383.
- Hannon E, et al. (2016). *Nature Neuroscience*, 19, 48-54.
- Kang HJ, et al. (2011). *Nature*, 478, 483-489.
- O'Brien HE, et al. (2018). *Genome Biology*, 19, 194.
- Polioudakis D, et al. (2019). *Neuron*, 103, 785-801.
- Trubetskoy V, et al. (2022). *Nature*, 604, 502-508.
- Wallace C (2021). *PLoS Genetics*, 17, e1009440.
- Zhu Z, et al. (2016). *Nature Genetics*, 48, 481-487.

## Licence

MIT. See `LICENSE`.
