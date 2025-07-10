# This script is used to calculate the alternate allele frequencies of specific positions for each population in 1000 Genomes Project
# This requires a prepared positions.txt file as an input. This file should have 4 columns (Chr, pos, REF, ALT). 
# I used tabs to separate columns, and I did not use headers. 
# This script also requires an integrated_call_samples_v3.20130502.ALL.panel from the 1000 Genomes Project. 
# The output will list all human populations in alphabetical order, chr, pos, REF, ALT, and ALT allele frequency for that population. 


#Extract your list of human samples
awk 'NR>1{print $2}' /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/integrated_call_samples_v3.20130502.ALL.panel | sort | uniq > populations.txt

while read pop; do
  awk -v pop=$pop 'NR>1 && $2==pop{print $1}' /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/integrated_call_samples_v3.20130502.ALL.panel > ${pop}.samples
done < populations.txt

#Calculate ALT frequencies for each population and position
while read pop; do
  sample_file="${pop}.samples"
  # Skip if sample list is empty
  if [ ! -s "$sample_file" ]; then continue; fi

  while read chr pos ref alt; do
    vcf_file="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    # Extract the variant for this position and population
    bcftools view -S $sample_file -r ${chr}:${pos}-${pos} $vcf_file | \
      bcftools query -f "${pop}\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n"
  done < positions.txt
done < populations.txt > all_populations_allele_freqs.txt
