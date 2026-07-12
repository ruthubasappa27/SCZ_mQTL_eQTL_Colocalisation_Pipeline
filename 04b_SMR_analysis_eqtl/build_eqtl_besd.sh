#!/usr/bin/env bash
# Build per-gene BESD files for eQTL SMR analysis
# Custom BESD construction required (no pre-built BESD available for O'Brien eQTL,
# unlike FB_Brain_2 for mQTL). ESI files must have exactly 6 columns
# (CHR SNP CM POS A1 A2) -- an extra 7th MAF column causes --make-besd to fail silently.
# Author: Ruthushree P Basappa
# Date: June 2026
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
for gene in CYP2D6 NAGA SMDT1 WBP2NL HSPA9 KRT18P46; do
    echo "Building BESD for ${gene}..."
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
        "${project_root}/references/smr_eqtl/${gene}_harm.esi" > \
        "${project_root}/references/smr_eqtl/${gene}_harm_fix.esi"
    mv "${project_root}/references/smr_eqtl/${gene}_harm_fix.esi" \
       "${project_root}/references/smr_eqtl/${gene}_harm.esi"
    "${project_root}/references/smr/smr_x86" \
        --eqtl-summary "${project_root}/references/smr_eqtl/${gene}_harm.txt" \
        --make-besd \
        --out "${project_root}/references/smr_eqtl/${gene}_harm"
done
echo "All per-gene BESD files built"
