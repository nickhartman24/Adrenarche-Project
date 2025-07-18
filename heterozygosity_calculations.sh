"""
This script has been created to quickly calculate heterozygosity rates for individuals in the 1000 Genomes Project for a select set of genes.
Input required for this script is a VCF of modern human genomes in the 1000 Genomes Project. 
This script provides outputs at multiple steps, including:
- homozygosity rates for all individuals,
- heterozygosity rates for all individuals, 
- average homozygosity rate (one number), 
- heterozygosity rates averaged by superpopulation (AFR, AMR, EUR, EAS, SAS). 
"""

#!/bin/bash

# Format: "GENE:CHR:START-END"
GENES=(
  "CYB5A:18:71920819-71959110"
  "HSD3B1:1:120049833-120057677"
  "HSD3B2:1:119957773-119965657"
  "CYP17A1:10:104590288-104597170"
  "SULT2A1:19:48373724-48389572"
)

for gene_info in "${GENES[@]}"; do
  GENE=$(echo $gene_info | cut -d: -f1)
  CHR=$(echo $gene_info | cut -d: -f2)
  REGION=$(echo $gene_info | cut -d: -f3)
  
VCF="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
PANEL_FILE="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/integrated_call_samples_v3.20130502.ALL.panel"
VCF_OUT="${GENE}_1000GP_gene_region.vcf.gz"

# Get sample-to-superpopulation mapping once
awk '{print $1 "\t" $3}' "$PANEL_FILE" > sample_superpop.txt
sort sample_superpop.txt > sample_superpop.sorted.txt

for gene_info in "${GENES[@]}"; do
  GENE=$(echo $gene_info | cut -d: -f1)
  CHR=$(echo $gene_info | cut -d: -f2)
  REGION=$(echo $gene_info | cut -d: -f3)

  # 1. Extract gene region
  bcftools view -r "${CHR}:${REGION}" "$VCF" -Oz -o "$VCF_OUT"

  # 2. Calculate per-individual HOMOzygosity
  # Columns in .het are: INDV O(HOM) E(HOM) N(NM) F
  vcftools --gzvcf "$VCF_OUT" --het --out "${GENE}_1000GP_homozyg_freq"

  # 3. Calculate individual heterozygosity
  awk 'NR>1 {print $1, ($4-$2)/$4}' "${GENE}_1000GP_homozyg_freq.het" > "${GENE}_1000GP_indv_het_freq"

  # 4. Calculate average heterozygosity (all samples)
  awk 'NR>1 {sum += ($4-$2)/$4; n++} END {print "Average heterozygosity:", sum/n}' "${GENE}_1000GP_homozyg_freq.het" > "${GENE}_average_heterozygosity_ALL.txt"

  # 5. Merge sample IDs with superpop
  sort "${GENE}_1000GP_indv_het_freq" > "${GENE}_1000GP_indv_het.sorted.txt"
  join -1 1 -2 1 "${GENE}_1000GP_indv_het.sorted.txt" sample_superpop.sorted.txt > "${GENE}_het_rates_superpop.txt"

  # 6. Average heterozygosity by superpopulation
  awk '{count[$3]++; sum[$3]+=$2} END {for (pop in sum) print pop, sum[pop]/count[pop]}' "${GENE}_het_rates_superpop.txt" > "${GENE}_heterozygosity_by_superpopulation.txt"
done

