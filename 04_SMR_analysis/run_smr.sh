#!/usr/bin/env bash
# Run SMR analysis for all 6 colocalised loci
# Author: Ruthushree P Basappa
# Date: May 2026
# Reference: Zhu et al. 2016, Nature Genetics
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
mkdir -p "${project_root}/results/smr"
"${project_root}/references/smr/smr_x86" \
  --bfile "${project_root}/references/g1000_eur_chrbp/g1000_eur" \
  --gwas-summary "${project_root}/data/gwas_smr_chrbp.txt" \
  --beqtl-summary "${project_root}/references/FB_Brain_2/FB_Brain_2" \
  --out "${project_root}/results/smr/smr_all_genes" \
  --extract-probe "${project_root}/04_SMR_analysis/all_probes.txt" \
  --peqtl-smr 5e-8 \
  --diff-freq 0.4 \
  --diff-freq-prop 0.6 \
  --thread-num 4
