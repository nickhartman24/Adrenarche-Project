#WHAT I'M DOING HERE - NUCLEOTIDES NOT FOUND IN AFR POPULATIONS (PRESUMED TO BE ANCESTRAL) BUT ARE PRESENT IN NEANDERTHALS. using slurm. example here is using SULT2A1, including UTRs. 

#starting with a slurm script


#!/bin/bash
# number of nodes
# we're doing one task per node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# Request to run on Villanea Lab node
#SBATCH --qos=blanca-villanea
##SBATCH --qos=preemptable
# Request runtime:
#SBATCH --time=01-00:00:00
# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
##SBATCH --mem=8G
# Specify a job name:
#SBATCH --job-name slurm_cyp17_bcftools
# Specify an output file
#SBATCH --output slurm_cyp17_bcftools.out
# Run a command
module load bcftools

## output files for each archaic, containing just the vcf data for the gene region
bcftools view -Oz /pl/active/villanea_lab/data/archaic_data/altai_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr19_mq25_mapab100.vcf.gz -r 19:48373724-48389572 -o sult2a_altai_Jan2025.vcf.gz
bcftools index sult2a_altai_Jan2025.vcf.gz

bcftools view -Oz /pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr19_mq25_mapab100.vcf.gz -r 19:48373724-48389572 -o sult2a_vindija_Jan2025.vcf.gz
bcftools index sult2a_vindija_Jan2025.vcf.gz

bcftools view -Oz /pl/active/villanea_lab/data/archaic_data/chagyrskaya_VCF/ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr19.noRB.vcf.gz -r 19:48373724-48389572 -o sult2a_chagyrskaya_Jan2025.vcf.gz
bcftools index sult2a_chagyrskaya_Jan2025.vcf.gz

## merge the archaic gene regions into one file
bcftools merge sult2a_altai_Jan2025.vcf.gz sult2a_vindija_Jan2025.vcf.gz sult2a_chagyrskaya_Jan2025.vcf.gz -Oz -o sult2a_merged_neanderthals_Jan2025.vcf.gz
bcftools index sult2a_merged_neanderthals_Jan2025.vcf.gz

## output a file for only YRI population for the gene region
bcftools view -Oz /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -r 19:48373724-48389572 -S YRI_samples_IDs.txt -o sult2a_YRI_Jan2025.vcf.gz 
bcftools index sult2a_YRI_Jan2025.vcf.gz

## Filter the merged archaic file for snps
bcftools filter -r 19:48373724-48389572 -i 'TYPE="snp"' sult2a_merged_neanderthals_Jan2025.vcf.gz -Oz -o sult2a_nea_snps_Jan2025.vcf.gz
bcftools index sult2a_nea_snps_Jan2025.vcf.gz

## removing headers and filtering to only archaic snps chromosome and position outputs a file that is used in the next step 
bcftools query -f '%CHROM\t%POS\n' sult2a_nea_snps_Jan2025.vcf.gz | awk '{print $1"\t"$2}' > sult2a_nea_snps_chr_pos_Jan2025.txt

## using the previously created file, filter the full YRI gene region down to only positions that align with the archaic snps
bcftools view -Oz -T sult2a_nea_snps_chr_pos_Jan2025.txt sult2a_YRI_Jan2025.vcf.gz -o filtered_sult2a_YRI_to_nea_snp_pos_Jan2025.vcf.gz
bcftools index filtered_sult2a_YRI_to_nea_snp_pos_Jan2025.vcf.gz

## merge the data for the archaics and the YRI individuals for the archaic snp positions
bcftools merge sult2a_nea_snps_Jan2025.vcf.gz filtered_sult2a_YRI_to_nea_snp_pos_Jan2025.vcf.gz -Oz -o sult2a_snps_merged_nea_YRI_Jan2025.vcf.gz
bcftools index sult2a_snps_merged_nea_YRI_Jan2025.vcf.gz

vcftools --gzvcf filtered_sult2a_YRI_to_nea_snp_pos_Jan2025.vcf.gz --freq --out sult2a_YRI_freq_at_nea_snp_pos

vcftools --gzvcf sult2a_nea_snps_Jan2025.vcf.gz --freq --out sult2a_freq_nea_snps
