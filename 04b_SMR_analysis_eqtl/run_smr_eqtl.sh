#!/usr/bin/env bash
# Run SMR analysis for eQTL colocalised gene signals
# Author: Ruthushree P Basappa
# Date: June 2026
# Reference: Zhu et al. 2016, Nature Genetics
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
mkdir -p "${project_root}/results_eqtl/smr"

"${project_root}/references/smr/smr_x86" \
  --bfile "${project_root}/references/g1000_eur_chrbp/g1000_eur" \
  --gwas-summary "${project_root}/data/gwas_smr_chrbp.txt" \
  --eqtl-flist "${script_dir}/eqtl_flist.txt" \
  --out "${project_root}/results_eqtl/smr/smr_all_genes_eqtl" \
  --peqtl-smr 1e-4 \
  --disable-freq-ck \
  --thread-num 4
