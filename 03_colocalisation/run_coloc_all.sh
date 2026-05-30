#!/bin/bash
cd ~/Desktop/Thesis

tail -n +2 COLOC-reporter/scripts/set.regions.txt | while IFS=$'\t' read traits region; do
    locus=$(echo $region | tr ',' '_')
    outdir="COLOC-reporter/coloc/results/${locus}"
    
    # Skip if already completed
    if [ -f "${outdir}_coloc/colocalisation.log" ]; then
        echo "Skipping locus $region - already completed"
        continue
    fi
    
    echo "Running locus: $region"
    Rscript COLOC-reporter/scripts/colocaliseRegion.R \
      --plink ~/Desktop/Thesis/plink_mac_20250819/plink \
      --traits $traits \
      --set_locus $region \
      --LDreference ~/Desktop/Thesis/COLOC-reporter/ld_reference/EUR_phase3_chr \
      --GWASconfig ~/Desktop/Thesis/COLOC-reporter/scripts/GWAS_samples.txt \
      --runMode trySusie \
      --out ~/Desktop/Thesis/COLOC-reporter/coloc/results/${locus} \
      --scriptsDir ~/Desktop/Thesis/COLOC-reporter/scripts
done

echo "ALL LOCI COMPLETED!"
