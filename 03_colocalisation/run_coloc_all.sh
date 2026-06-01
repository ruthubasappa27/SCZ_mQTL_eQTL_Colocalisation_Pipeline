#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

tail -n +2 "${PROJECT_ROOT}/03_colocalisation/set.regions.txt" | while IFS=$'\t' read -r traits region; do
    locus=$(echo "${region}" | tr ',' '_')
    outdir="${PROJECT_ROOT}/results/coloc/${locus}"

    if [ -f "${outdir}_coloc/colocalisation.log" ]; then
        echo "Skipping locus ${region} - already completed"
        continue
    fi

    echo "Running locus: ${region}"
    Rscript "${PROJECT_ROOT}/COLOC-reporter/scripts/colocaliseRegion.R" \
      --plink "${PROJECT_ROOT}/plink_mac_20250819/plink" \
      --traits "${traits}" \
      --set_locus "${region}" \
      --LDreference "${PROJECT_ROOT}/COLOC-reporter/ld_reference/EUR_phase3_chr" \
      --GWASconfig "${PROJECT_ROOT}/03_colocalisation/GWAS_samples.txt" \
      --runMode trySusie \
      --out "${PROJECT_ROOT}/results/coloc/${locus}" \
      --scriptsDir "${PROJECT_ROOT}/COLOC-reporter/scripts"
done

echo "ALL LOCI COMPLETED!"

