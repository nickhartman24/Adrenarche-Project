#DESCRIPTION OF WHAT I'M DOING HERE. did this in a mac terminal, but in a conda environment. used bash throughout. also include that it creates multiple different types of outputs, and link to og haplostrips. Explain what inputs I'm using (1000GP chr18 vcf, vindija vcf, 0.001 threshold, etc) 
#WHAT THIS REQUIRES BEFORE STARTING. THINK ABOUT: INSTALLING HAPLOSTRIPS, CONDA ENVIRONEMT, ETC. SAMPLES_HAPLO FILE AND IT'S LAYOUT  

#conduct haplostrips on gene CYB5A. Split into 3 groups because the figure creator aspect of Haplostrips is limited to a certain amount of populations allowed. 
#using the positions from UCSC genome browser, coding region only, no UTRs included. Nea (Vindija Neanderthal) must be included each round as the reference.
#This script is being used as a test for introgression, and African populations are not expected to be introgressed. The only 1000 Genomes African population used is YRI. ASW and ACB are also excluded becuase of their African Ancestry. 

python haplostrips.py -v /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -v /pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr18_mq25_mapab100.vcf.gz -P samples_haplo -i 18:71920819-71959110 -p Nea,CEU,TSI,GBR,GIH,PJL,BEB -c 0.001 -a -t -o cyb5a_first_group

python haplostrips.py -v /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -v /pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr18_mq25_mapab100.vcf.gz -P samples_haplo -i 18:71920819-71959110 -p Nea,IBS,ITU,STU,FIN,PUR,CLM -c 0.001 -a -t -o cyb5a_second_group

python haplostrips.py -v /pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -v /pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr18_mq25_mapab100.vcf.gz -P samples_haplo -i 18:71920819-71959110 -p Nea,MXL,PEL,CDX,KHV,CHB,CHS,JPT,YRI -c 0.001 -a -t -o cyb5a_third_group

#combine the four group files into one. 
head -n 1 $(ls cyb5a*group.distances_tab | head -n 1) > cyb5a_merged.distances_tab
tail -n +2 -q cyb5a*group.distances_tab >> cyb5a_merged.distances_tab

#Sort by the third column, which is distance, while preserving the header.  
(head -n 1 cyb5a_merged.distances_tab && tail -n +2 cyb5a_merged.distances_tab | sort -k3,3n) > cyb5a_sorted.distances_tab
