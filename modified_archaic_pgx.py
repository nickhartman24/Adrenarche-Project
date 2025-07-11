"""
The original version of this script was written by Kelsey Witt in 2021 (https://github.com/kelsey-witt/archaic_pgx/tree/main). 
This script contains modified versions of https://github.com/kelsey-witt/archaic_pgx/blob/main/archaic_allele_freqs_pharmacogenes_mp.py and 
https://github.com/kelsey-witt/archaic_pgx/blob/main/compare_rare_afr_archaic.py. 

The first part of this script is titled "##Archaic Allele Frequencies", and is a modified 'archaic_allele_freqs_pharmacogenes_mp.py' script. 
This part of this script has two purposes. The first uses the specified adrenarche genes and their coordinates
to identify single nucleotide variants for each gene that are at low frequency in Africa, and the output of this is Gene + "_AFRrare_snp_freqs.csv".
The second purpose is to output a file with all genotypes where at least one archaic individual has an alternate allele. This output is 
"archaic_SNPs_at_rare_pos_" + Gene + ".csv"

The last portion of this script is the modified "compare_rare_afr_archaic.py", under the section "##Compare_rare_afr_archaic". It takes the two 
outputs from the first part of the script, and combines them to identify SNVs that are rare in Africa and shared with archaic humans. The output of this is 
afile called "archaic_mp_nonafr_snp_freqs.csv". 
This is a way to determine at what positions Neanderthals and the ancestral human allele are different, and could be candidates for alleles 
with functional differences.
"""

#Inputs requires archaic subsets for each gene. Did Vindija, Altai, Chagyrskaya, and Denisovan. 
# This is the example for HSD3B1, and need to be repeated for each gene of interest. 
# After this, the script is automated. 
bcftools view -r 1:120049833-120057677 -o Altai_hsd3b1.vcf.gz -O z /pl/active/villanea_lab/data/archaic_data/altai_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr1_mq25_mapab100.vcf.gz

bcftools view -r 1:120049833-120057677 -o Vindija_hsd3b1.vcf.gz -O z /pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr1_mq25_mapab100.vcf.gz

bcftools view -r 1:120049833-120057677 -o Chagyrskaya_hsd3b1.vcf.gz -O z /pl/active/villanea_lab/data/archaic_data/chagyrskaya_VCF/ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr1.noRB.vcf.gz

bcftools view -r 1:120049833-120057677 -o Denisovan_hsd3b1.vcf.gz -O z /pl/active/villanea_lab/data/archaic_data/denisova_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr1_mq25_mapab100.vcf.gz

##Archaic Allele Frequencies 
import gzip
import re

adrenGenes = {"hsd3b1": ["1",120049833,120057677], "hsd3b2": ["1",119957773,119965657], "cyp17a1": ["10",104590288,104597170], "sult2a1": ["19",48373724,48389572], "cyb5a": ["18",71920819,71959110]}

pop_file = "/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/integrated_call_samples_v3.20130502.ALL.panel"

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
afrPops = ["YRI", "LWK", "GWD", "MSL", "ESN"]
archaics = ["Denisovan","Altai","Chagyrskaya","Vindija"]
Pop_Tracker = {}
rareAfrAlleles = []
archaicGenotypes = {}

for pop in populations:
    Pop_Tracker[pop] = []

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] == "AFR" and line_col[1] in afrPops: 
            Pop_Tracker["AFR"].append(col_counter)
        elif line_col[2] != "AFR":
            pop = line_col[1]
            Pop_Tracker[pop].append(col_counter)
        col_counter += 1

def decode_split_line(line):
    #line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

def snp_is_biallelic(line):
    ref_allele = line[3]
    alt_allele = line[4]
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return True
    else:
        return False

def parse_arch_line(line):
    position, _, ref_allele, arch_allele = line[1:5]
    if arch_allele == '.':
        arch_allele = ref_allele
    arch_info = line[9]
    return [position, ref_allele, arch_allele, arch_info]

def parse_mod_line(line):
    position, rsVal, ref_allele, mod_alt_allele = line[1:5]
    mod_info = line[7]
    return [position, ref_allele, mod_alt_allele, mod_info, rsVal]

def calc_afr_allele_freq(line):
    afr_allele_total = len(Pop_Tracker["AFR"]) * 2
    afr_allele_1_sum = 0
    for afr_ind in Pop_Tracker["AFR"]:
        afr_allele_1_sum += line[afr_ind].count("1")
    afr_allele_1_fq = float(afr_allele_1_sum/afr_allele_total)
    afr_allele_1_round = round(afr_allele_1_fq)
    return afr_allele_1_round

def identify_non_afr_alleles(line):
    non_afr_allele = "null"
    afr_alt_af = calc_afr_allele_freq(line)
    if afr_alt_af < 0.01 or afr_alt_af > 0.99:
        if afr_alt_af < 0.01: #allele of interest is alternate
            non_afr_allele = line[4]
        else: #allele of interest is reference
            non_afr_allele = line[3]
    return non_afr_allele

def arch_passes_quality_check(arch_info):
    arch_data_parse = arch_info.split(":")
    if len(arch_data_parse) == 11: #Altai Nea or Den info line
        genotype_quality = arch_data_parse[2]
        mapping_quality = arch_data_parse[9]
        if genotype_quality == ".":
            genotype_quality = 0.0
        if mapping_quality == ".":
            mapping_quality = 0.0
        if float(genotype_quality) >= 40.0 and float(mapping_quality) >= 30.0:
            return True
        else:
            return False
    elif len(arch_data_parse) == 8: #Vindija Nea info line
        genotype_quality = arch_data_parse[7]
        if float(genotype_quality) >= 40.0:
            return True
        else:
            return False

def calc_pop_freq(population, nonafr_allele):
    pop_allele_count = 0
    allele_total = len(Pop_Tracker[population]) * 2
    for ind in Pop_Tracker[population]:
        pop_allele_count += mod_line[ind].count(nonafr_allele)
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return rounded_freq

for adrenGene in adrenGenes:
    rareAfrAlleles = []
    archaicGenotypes = {}
    outfile = adrenGene + "_AFRrare_snp_freqs.csv"
    f = open(outfile, 'w')
    headercols = ["adrenacogene","chromosome","position","rs value","ref","alt","rare allele"] + populations
    headerline = ",".join(headercols)+"\n"
    f.write(headerline)
    chromosome, adrenSt, adrenEnd = adrenGenes[adrenGene][0:3]
    modern_file = "/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    with gzip.open(modern_file, 'rt') as mf:
        for mod_line in mf:
            mod_line = decode_split_line(mod_line)
            #print(mod_line)
            if '#' not in mod_line:
                mod_vars = parse_mod_line(mod_line)
                mod_position, mod_ref, mod_alt = mod_vars[0:3]
                mod_rs = mod_vars[4]
                if int(mod_position) >= adrenSt and int(mod_position) <= adrenEnd:
                    freq_by_pop = []
                    if snp_is_biallelic(mod_line):
                        mod_alt_allele = mod_line[4]
                        for pop in populations[1:]:
                                alleleCounter = 0
                                allele_total = len(Pop_Tracker[pop]) * 2
                                for ind in Pop_Tracker[pop]:
                                    alleleCounter += mod_line[ind].count("1")
                        non_afr_allele = identify_non_afr_alleles(mod_line)
                        if non_afr_allele != "null":
                            if non_afr_allele == mod_ref:
                                non_afr_genotype = "0"
                            elif non_afr_allele == mod_alt:
                                non_afr_genotype = "1"
                            rareAfrAlleles.append(int(mod_position))
                            archaicGenotypes[int(mod_position)] = ["*","*","*","*",mod_ref,mod_alt]
                            for pop in populations:
                                arch_freq = calc_pop_freq(pop, non_afr_genotype)
                                freq_by_pop.append(str(arch_freq))
                            datacols = [adrenGene,chromosome,mod_position,mod_rs,mod_ref,mod_alt,non_afr_genotype]+freq_by_pop
                            dataline = ",".join(datacols) +"\n"
                            f.write(dataline)
    f.close()
    for i in range(0,len(archaics)):
        archaic = archaics[i]
        archInfile = archaic + "_" + adrenGene + ".vcf.gz"
        archOutfile = "archaic_SNPs_at_rare_pos_" + adrenGene + ".csv"
        f = open(archOutfile, 'w')
        headercols = ["adrenacogene","chromosome","position","Den","Altai","Chag","Vind","ref","alts"]
        headerline = ",".join(headercols)+"\n"
        f.write(headerline)
        with gzip.open(archInfile, 'rt') as af:
            for archline in af:
                if '#' not in archline:
                    archline = decode_split_line(archline)
                    arch_vars = parse_arch_line(archline) 
                    arch_position, ref_allele, arch_alt_allele, arch_info = arch_vars[0:4]
                    if int(arch_position) in rareAfrAlleles:
                        if arch_alt_allele not in archaicGenotypes[int(arch_position)][5]:
                            archaicGenotypes[int(arch_position)][5] += "_" + arch_alt_allele
                        arch_genotype = arch_info[0:3]
                        archaicGenotypes[int(arch_position)][i] = arch_genotype
    for position in rareAfrAlleles:
        pos = int(position)
        datacols = [adrenGene,chromosome,str(pos)]+archaicGenotypes[pos]
        dataline = ",".join(datacols) +"\n"
        f.write(dataline)
    f.close()


##Compare_rare_afr_archaic
adrenGenes = {"hsd3b1”,“hsd3b2","cyp17a1","sult2a1","cyb5a"}

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]

outfile = "archaic_mp_nonafr_snp_freqs.csv"
f = open(outfile, 'w')
headercols = ["adren_gene","chromosome","position","rs value","ref","alt","rare allele",”Den”,"Altai","Chagyrskaya","Vindija","arch_alt"] + populations
headerline = ",".join(headercols)+"\n"
f.write(headerline)
for adrenGene in adrenGenes:
    afrInfile = adrenGene + "_AFRrare_snp_freqs.csv"
    print(afrInfile)
    archInfile = "archaic_SNPs_at_rare_pos_" + adrenGene + ".csv"
    h = open(archInfile)
    with open(afrInfile) as g:
        h.readline()
        next(g)
        for afrLine in g:
            afrSharedArchaic = False
            archLine = h.readline()
            archSpline = archLine.split(",")
            denGeno,altGeno,chagGeno,vinGeno = archSpline[3:7]
            if denGeno != "*" or altGeno != "*" or chagGeno != "*" or vinGeno != "*": #ensures there's data
                archPos = archSpline[2]
                afrSpline = afrLine.split(sep=",")
                afrPos = afrSpline[2]
                if archPos == afrPos:
                    nonAfrGeno = afrSpline[6]
                    for archaic in (denGeno,altGeno,chagGeno,vinGeno):
                        if nonAfrGeno in archaic:
                            afrSharedArchaic = True
                    if afrSharedArchaic:
                        afrStart=afrSpline[0:7]
                        afrEnd=afrSpline[7:]
                        archAlt=archSpline[8][:-1]
                        outLine = afrStart+[denGeno,altGeno,chagGeno,vinGeno,archAlt]+afrEnd
                        dataline = ",".join(outLine) +"\n"
                        f.write(dataline)   
f.close()
