# This is a slurm script that automates the process of identifying Neanderthal alleles for specified genes. 
# We use Yoruban (YRI) genomes from the 1000 Genomes Project, making the assumption that YRI has no archaic introgression, 
# and that the REF allele is ancestral to modern humans if all YRI share the REF allele. 
# We compare this to the three Neanderthal genomes. Alleles are determined to be "Unique Neanderthal variants" if all Neanderthals share the ALT allele 
#(ALT freq=1), at positions where all YRI share the REF allele (REF freq=1). 
# Inputs include VCFs for the three Neanderthal genomes, the integrated VCF with all 1000 Genomes Project modern human samples, and a list of the YRI sample IDs. 
# Outputs include YRI .frq files and Neanderthal .frq files for snps in each specified gene. 


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=blanca-villanea
#SBATCH --time=01-00:00:00
#SBATCH --mem=8G
#SBATCH --job-name slurm_nea_african_bcftools
#SBATCH --output slurm_nea_african_bcftools_%j.out

module load bcftools
module load vcftools

# Format: "GENE:CHR:START-END"
REGIONS=(
  "HSD3B1:1:120049833-120057677"
  "HSD3B2:1:119957773-119965657"
  "CYB5A:18:71920819-71959110"
  "CYP17A1:10:104590288-104597170"
  "SULT2A1:19:48373724-48389572"
)

for region in "${REGIONS[@]}"; do
  GENE=$(echo $region | cut -d: -f1)
  CHR=$(echo $region | cut -d: -f2)
  COORDS=$(echo $region | cut -d: -f3)

# Inputs: VCF of each Neanderthal genome, and the modern human 1000 Genomes Project integrated genomes. A text file of the Yoruban (YRI) sample ids
ALTAI_VCF="/pl/active/villanea_lab/data/archaic_data/altai_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${CHR}_mq25_mapab100.vcf.gz"
VINDIJA_VCF="/pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${CHR}_mq25_mapab100.vcf.gz"
CHAGYRSKAYA_VCF="/pl/active/villanea_lab/data/archaic_data/chagyrskaya_VCF/ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${CHR}.noRB.vcf.gz"
MODERN_VCF="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
YRI_LIST="YRI_samples_IDs.txt"

 # 1. Extract archaic vcf lines for each gene, and index
  bcftools view -Oz "$ALTAI_VCF" -r ${CHR}:${COORDS} -o ${GENE}_altai.vcf.gz
  bcftools index ${GENE}_altai.vcf.gz

  bcftools view -Oz "$VINDIJA_VCF" -r ${CHR}:${COORDS} -o ${GENE}_vindija.vcf.gz
  bcftools index ${GENE}_vindija.vcf.gz

  bcftools view -Oz "$CHAGYRSKAYA_VCF" -r ${CHR}:${COORDS} -o ${GENE}_chagyrskaya.vcf.gz
  bcftools index ${GENE}_chagyrskaya.vcf.gz

  # 2. Merge archaic regions
  bcftools merge ${GENE}_altai.vcf.gz ${GENE}_vindija.vcf.gz ${GENE}_chagyrskaya.vcf.gz -Oz -o ${GENE}_merged_neanderthals.vcf.gz
  bcftools index ${GENE}_merged_neanderthals.vcf.gz

  # 3. From the full modern human vcf, filter to just YRI samples, and index
  bcftools view -Oz "$MODERN_VCF" -r ${CHR}:${COORDS} -S "$YRI_LIST" -o ${GENE}_YRI.vcf.gz
  bcftools index ${GENE}_YRI.vcf.gz

  # 4. Filter merged archaic for SNPs
  bcftools filter -r ${CHR}:${COORDS} -i 'TYPE="snp"' ${GENE}_merged_neanderthals.vcf.gz -Oz -o ${GENE}_nea_snps.vcf.gz
  bcftools index ${GENE}_nea_snps.vcf.gz

  # 5. Output chr/pos list for archaic SNPs, will be used in next step
  bcftools query -f '%CHROM\t%POS\n' ${GENE}_nea_snps.vcf.gz | awk '{print $1"\t"$2}' > ${GENE}_nea_snps_chr_pos.txt

  # 6. Filter YRI VCF to only archaic SNP positions
  bcftools view -Oz -T ${GENE}_nea_snps_chr_pos.txt ${GENE}_YRI.vcf.gz -o filtered_${GENE}_YRI_to_nea_snp_pos.vcf.gz
  bcftools index filtered_${GENE}_YRI_to_nea_snp_pos.vcf.gz

  # 7. Merge archaic SNPs with filtered YRI. Not used for allele frequency, but is generally good to have this file with full vcf lines. 
  bcftools merge ${GENE}_nea_snps.vcf.gz filtered_${GENE}_YRI_to_nea_snp_pos.vcf.gz -Oz -o ${GENE}_snps_merged_nea_YRI.vcf.gz
  bcftools index ${GENE}_snps_merged_nea_YRI.vcf.gz

  # 8. Allele frequencies for YRI and Neanderthals. "Unique Neanderthal variants" will be positions where YRI have a REF freq=1, and Nea have ALT freq=1. 
  vcftools --gzvcf filtered_${GENE}_YRI_to_nea_snp_pos.vcf.gz --freq --out ${GENE}_YRI_freq_at_nea_snp_pos
  vcftools --gzvcf ${GENE}_nea_snps.vcf.gz --freq --out ${GENE}_freq_nea_snps

done
