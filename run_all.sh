#!/bin/bash
# Runs the whole pipeline in order. Usage: bash run_all.sh
set -e
mkdir -p logs
Rscript scripts/00_check_setup.R 2>&1 | tee logs/00_check_setup.txt
Rscript scripts/01_qc_gwas.R
Rscript scripts/02_define_loci.R
Rscript scripts/03_qc_mqtl.R
Rscript scripts/04_qc_eqtl.R
Rscript scripts/05_colocalisation.R mqtl
Rscript scripts/06_results.R mqtl
Rscript scripts/05_colocalisation.R eqtl
Rscript scripts/06_results.R eqtl
Rscript scripts/07_smr.R mqtl
Rscript scripts/07_smr.R eqtl
Rscript scripts/08_annotate_genes.R
Rscript scripts/09_audit_direction.R
Rscript scripts/10_verify_rsids.R
Rscript scripts/11_prioritise.R
Rscript scripts/12_developmental_expression.R
Rscript scripts/13_cell_type_expression.R
Rscript scripts/14_main_tables.R
Rscript scripts/15_workflow_diagram.R
Rscript scripts/16_flow_diagram.R
Rscript scripts/17_main_figures.R
Rscript scripts/18_supplementary_workbook.R
