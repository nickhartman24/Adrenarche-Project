# This script uses HAPLOSTRIPS (https://bitbucket.org/dmarnetto/haplostrips/src/master/) 
# as a preliminary test for Neanderthal introgression in a select group of genes. 
# Written and executable in bash/shell. I conducted this in a conda environment
# with Haplostrips software package (version 1.3; Marnetto & Huerta-Sánchez, 2017) implemented in Python 2.7.18. 
# Haplostrips will produce many output files (.haps, .pdf, .dist.pdf, .distances_tab, .mat). 
# The .pdf file creates a figure visualizing haplotypes of populations compared to a reference. 
# necessary inputs include VCFs for each population. I used an integrated VCF for modern humans from the 1000 Genome Project, 
# and a VCF for the Vindija Neanderthal. This script also requires a "samples_id.txt" (I called "samples_haplo") 
# with all individuals(humans and Neanderthals) that you want to analyze. 
# Samples_id.txt has two columns (sample_id and pop_code), tab separated and including headers. 
# The .pdf creating aspect of Haplostrips is limited so a certain volume of populations, which is why this script splits the human populations into 
# 3 groups, with Vindija Neanderthal in each run as the reference genome. 
# Because African populations are not expected to be introgressed, I only used Yorubans (YRI) to reduce noise.
# I also did not include ASW and ACB populations for the same reason. 
# The last portion of this script combines the .distances_tab files for the 3 groups into one, and sorts individuals by distance from Neanderthal.  


#!/bin/bash

# Define your list of genes: format is "GENE_NAME:CHR:START-END"
GENES=(
  "HSD3B1:1:120049833-120057677"
  "HSD3B2:1:119957773-119965657"
  "CYB5A:18:71920819-71959110"
  "CYP17A1:10:104590288-104597170"
  "SULT2A1:19:48373724-48389572"
)

for gene_info in "${GENES[@]}"; do
  GENE=$(echo $gene_info | cut -d: -f1)
  CHR=$(echo $gene_info | cut -d: -f2)
  REGION=$(echo $gene_info | cut -d: -f3)
  
# Population groups (customize as needed)
GROUP1="Nea,CEU,TSI,GBR,GIH,PJL,BEB"
GROUP2="Nea,IBS,ITU,STU,FIN,PUR,CLM"
GROUP3="Nea,MXL,PEL,CDX,KHV,CHB,CHS,JPT,YRI"

#VCF inputs
VCF_1000G="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
VCF_NEAN="/pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${CHR}_mq25_mapab100.vcf.gz"
SAMPLES="samples_haplo"
CUTOFF=0.001

  # Run Haplostrips for each group
  python haplostrips.py -v "$VCF_1000G" -v "$VCF_NEAN" -P "$SAMPLES" -i "${CHR}:${REGION}" -p "$GROUP1" -c "$CUTOFF" -a -t -o "${GENE}_first_group"
  python haplostrips.py -v "$VCF_1000G" -v "$VCF_NEAN" -P "$SAMPLES" -i "${CHR}:${REGION}" -p "$GROUP2" -c "$CUTOFF" -a -t -o "${GENE}_second_group"
  python haplostrips.py -v "$VCF_1000G" -v "$VCF_NEAN" -P "$SAMPLES" -i "${CHR}:${REGION}" -p "$GROUP3" -c "$CUTOFF" -a -t -o "${GENE}_third_group"

  # Merge output
  head -n 1 $(ls ${GENE}*group.distances_tab | head -n 1) > "${GENE}_merged.distances_tab"
  tail -n +2 -q ${GENE}*group.distances_tab >> "${GENE}_merged.distances_tab"

  # Sort merged output by distance
  (head -n 1 "${GENE}_merged.distances_tab" && tail -n +2 "${GENE}_merged.distances_tab" | sort -k3,3n) > "${GENE}_sorted.distances_tab"
done
