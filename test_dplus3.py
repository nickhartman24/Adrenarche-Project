import numpy as np
import gzip
import tabix
import argparse
import random
import re
import os
from collections import defaultdict

#FUNCTIONS
def get_names_of_populations_dict(string,populations=["pop1","pop2","pop3","outgroup"]):
  dict={}
  for pop in populations:
    search_string=r'{}\[(.*?)\]'.format(pop)
    population_names=re.findall(search_string,string)[0].split()
    dict[pop]=population_names #[header.index(ind) for ind in population_names]
  return(dict)

def get_index_populations_dict(vcf_name,string,populations=["pop1","pop2","pop3"]):
  #Get header of names
  with gzip.open(vcf_name,'rt') as file:
    for line in file:
      if line.startswith("#CHROM"):
        header=line.split()
        break
  #Get names in the population
  name_pop_dict=get_names_of_populations_dict(string,populations=populations)
  index_pop_dict=defaultdict()
  for key,values in name_pop_dict.items():
    index_pop_dict[key]=[header.index(ind) for ind in values]
  return index_pop_dict

def get_first_last_position_per_chromosome(infile):
    """Get first and last position in chromosome from dictionary"""
    first_last_position_chromosome_dict = {}
    with open(infile, "r") as file:
        for line in file:
            spline = line.strip().split(",")
            try:
                chrom = int(spline[0])      # Store as int for robust lookup
                first = int(spline[1])
                last = int(spline[2])
                first_last_position_chromosome_dict[chrom] = (first, last)
            except ValueError:
                continue
    return first_last_position_chromosome_dict

def get_ancestral_allele(pos,information_string):
  try:
    ancestral_allele=re.findall(r'AA=(.*?)\|;',information_string)[0].split("|")[0]
    ancestral_allele=ancestral_allele.upper()
  except IndexError:
    print(pos,information_string)
  return ancestral_allele

def get_allele_from_individual(reference_allele,alternative_allele,allele_integer):
  """Replace an individual allele '0' or '1' with the reference or alternative allele ['A','T','C','G']""" 
  allele_integer=str(allele_integer).strip()
  allele = allele_integer.replace('0',reference_allele).replace('1',alternative_allele)
  assert((allele==reference_allele) or (allele==alternative_allele)),"Sampled allele from an individual '{}' is not '0' or '1'".format(allele_integer)
  return allele

def allele_site_pattern(pop1_allele,pop2_allele,pop3_allele,pop4_allele):
  #BAAA or ABAA
  if(pop4_allele==pop3_allele):
    #BAAA
    if((pop3_allele==pop2_allele) & (pop1_allele!=pop2_allele)):
      return "BAAA"
    #ABAA
    elif((pop3_allele==pop1_allele) & (pop1_allele!=pop2_allele)):
      return "ABAA"
  #ABBA or BABA
  elif(pop4_allele!=pop3_allele):
    #ABBA
    if((pop3_allele==pop2_allele) & (pop1_allele!=pop2_allele)):
      return "ABBA"
    #ABAA
    elif((pop3_allele==pop1_allele) & (pop1_allele!=pop2_allele)):
      return "BABA"
  return None

def get_derived_allele(ancestral_allele,reference_allele,alternative_allele):
  if ancestral_allele==reference_allele:
    derived_allele=0
  elif ancestral_allele==alternative_allele:
    derived_allele=1
  return derived_allele

def get_derived_freq_per_population(array,ancestral_allele,reference_allele,alternative_allele):
  """Get derived allele frequency for each population. Input a dictionary with the alleles ["A","T","C","G"] for the biallelic
    site in each of the four populations. Output a dictionary with the float derived allele frequency for each population.
  """
  derived_allele=get_derived_allele(ancestral_allele=ancestral_allele,reference_allele=reference_allele,alternative_allele=alternative_allele)
  freq=(array.count(derived_allele),np.shape(array)[0])
  return freq

def update_introgression_stats_sums(p_one,p_two,p_three,p_four,list):
  dict=defaultdict(int)
  if "D" in list:
    #D statistic
    dict["d_numerator"]=((1-p_one)*p_two*p_three*(1-p_four))-(p_one*(1-p_two)*p_three*(1-p_four))
    dict["d_denominator"]=((1-p_one)*p_two*p_three*(1-p_four))+(p_one*(1-p_two)*p_three*(1-p_four))
  if "D+" in list:
    #D+ statistic
    dict["dplus_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_numerator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))-((1-p_one)*p_two*(1-p_three)*(1-p_four))
    dict["dplus_denominator"]=((1-p_one)*p_two*p_three*(1-p_four))+(p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_denominator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))+((1-p_one)*p_two*(1-p_three)*(1-p_four))
  if "Dancestral" in list:
    #D+ statistic
    dict["dancestral_numerator"]=(p_one*(1-p_two)*(1-p_three)*(1-p_four))-((1-p_one)*p_two*(1-p_three)*(1-p_four))
    dict["dancestral_denominator"]=(p_one*(1-p_two)*(1-p_three)*(1-p_four))+((1-p_one)*p_two*(1-p_three)*(1-p_four))
  if "fD" in list:
    #fD statistic
    dict["fd_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fD uses the max(p2,p3) as donor population
    pd=max(p_two,p_three)
    dict["fd_denominator"]=((1-p_one)*pd*pd*(1-p_four))-(p_one*(1-pd)*pd*(1-p_four))
  if "fDM" in list:
    #fDM statistic
    dict["fdm_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fDM has the same denominator as fD when p_two is greater than or equal to p3
  if p_one>=p_two:
    pd_m=max(p_two,p_three)
    dict["fdm_denominator"]=((1-p_one)*pd_m*pd_m*(1-p_four))-(p_one*(1-pd_m)*pd_m*(1-p_four))
  elif p_one<p_two:
    pd_m=max(p_one,p_three)
    dict["fdm_denominator"]=(pd_m*(1-p_two)*pd_m*(1-p_four))-((1-pd_m)*p_two*pd_m*(1-p_four))
  if "df" in list:
    #Df statsitic
    dict["df_numerator"]=((1-p_one)*p_two*p_three)-(p_one*(1-p_two)*p_three)
    dict["df_denominator"]=(2*p_one*p_two*p_three)+((1-p_one)*p_two*p_three)+(p_one*(1-p_two)*p_three)
  return(dict)

def calculate_introgression_stats(dict,statistics_list):
  results=defaultdict(int)
  if "D" in statistics_list:
    #D statistic
    try:
      results["D"]=(dict["ABBA"]-dict["BABA"])/(dict["ABBA"]+dict["BABA"])
    except ZeroDivisionError:
      results["D"]=np.nan
  if "D+" in statistics_list:
    #D+ statistic
    try:
      results["D+"]=((dict["ABBA"]-dict["BABA"])+(dict["BAAA"]-dict["ABAA"]))/(dict["ABBA"]+dict["BABA"]+dict["BAAA"]+dict["ABAA"])
    except ZeroDivisionError:
      results["D+"]=np.nan
  if "Dancestral" in statistics_list:
    #Dancestral statistic
    try:
      results["Dancestral"]=dict["dancestral_numerator"]/dict["dancestral_denominator"]
    except ZeroDivisionError:
      results["Dancestral"]=np.nan
  if "fD" in statistics_list:
    #fD statistic
    try:
      results["fD"]=dict["fd_numerator"]/dict["fd_denominator"]
    except ZeroDivisionError:
      results["fD"]=np.nan
  if "fDM" in statistics_list:
    #fDM statistic
    try:
      results["fDM"]=dict["fdm_numerator"]/dict["fdm_denominator"]
    except ZeroDivisionError:
      results["fDM"]=np.nan
  if "df" in statistics_list:
    #Df statsitic
    try:
      results["df"]=dict["df_numerator"]/dict["df_denominator"]
    except ZeroDivisionError:
      results["df"]=np.nan
    return results

  if "D+" in list:
    #D+ statistic
    dict["dplus_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_numerator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))-((1-p_one)*p_two*(1-p_three)*(1-p_four))
    dict["dplus_denominator"]=((1-p_one)*p_two*p_three*(1-p_four))+(p_one*(1-p_two)*p_three*(1-p_four))
    dict["dplus_denominator"]+=(p_one*(1-p_two)*(1-p_three)*(1-p_four))+((1-p_one)*p_two*(1-p_three)*(1-p_four))
  if "fD" in list:
    #fD statistic
    dict["fd_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fD uses the max(p2,p3) as donor population
    pd=max(p_two,p_three)
    dict["fd_denominator"]=((1-p_one)*pd*pd*(1-p_four))-(p_one*(1-pd)*pd*(1-p_four))
  if "fDM" in list:
    #fDM statistic
    dict["fdm_numerator"]=((1-p_one)*p_two*p_three*(1-p_four)) - (p_one*(1-p_two)*p_three*(1-p_four))
    #fDM has the same denominator as fD when p_two is greater than or equal to p3
    if p_one>=p_two:
      pd_m=max(p_two,p_three)
      dict["fdm_denominator"]=((1-p_one)*pd_m*pd_m*(1-p_four))-(p_one*(1-pd_m)*pd_m*(1-p_four))
    elif p_one<p_two:
      pd_m=max(p_one,p_three)
      dict["fdm_denominator"]=(pd_m*(1-p_two)*pd_m*(1-p_four))-((1-pd_m)*p_two*pd_m*(1-p_four))
  if "df" in list:
    #Df statsitic
    dict["df_numerator"]=((1-p_one)*p_two*p_three)-(p_one*(1-p_two)*p_three)
    dict["df_denominator"]=(2*p_one*p_two*p_three)+((1-p_one)*p_two*p_three)+(p_one*(1-p_two)*p_three)
  return(dict)

def get_introgressed_region_for_haplotype(file_name,individual):
  introgressed_tracts=[]
  with gzip.open(file_name,'rt') as file:
    for line in file:
      if line.startswith("##"):
        continue
      spline=line.split()
      if int(spline[1]) in individual:
        introgressed_tracts.append((int(spline[2]),int(spline[3])))
  return introgressed_tracts

def get_introgressed_region_for_frequencies(file_name, chromosome):
    introgressed_tracts = []
    with open(file_name, "r") as file:
        for line in file:
            spline = line.strip().split(",")
            if int(spline[0]) == chromosome:
                introgressed_tracts.append((int(spline[1]), int(spline[2])))
    return introgressed_tracts

def number_intro_bases(window_start,window_stop,intro_start,intro_stop):
  bases=0
  #Introgressed region 100% overlap
  if(window_start >= intro_start and window_stop <= intro_stop):
    bases=window_stop - window_start
  #Introgressed region overlap at beginning
  if(window_start >= intro_start and window_stop >= intro_stop and intro_stop > window_start):
    bases = intro_stop-window_start + 1
  #Introgressed region overlap at end
  if(window_start <= intro_start and window_stop <= intro_stop and intro_start < window_stop):
    bases = window_stop - intro_start + 1
  #Introgressed region within window
  if(window_start <= intro_start and window_stop >= intro_stop):
    bases = intro_stop - intro_start + 1
  return(bases)

#Parameters
parser=argparse.ArgumentParser()

parser.add_argument("-o","--outfile",type=str,help="Path to outfile directory.")
parser.add_argument("-human_infile",type=str,dest="human_vcf",help="Path to vcfs.")
parser.add_argument("-archaic_infile",type=str,dest="archaic_vcf",help="Path to vcfs.")
parser.add_argument("-introgressed_infile",type=str,help="Path to introgression map file.")
parser.add_argument("-chromosome",type=str,help="Chromosome number to analyze.")
parser.add_argument("--report",nargs='?',const=100,type=int,help="When to report progress on windows.")
parser.add_argument("-ws","--window_size",type=int,help="Size of the sliding non-overlapping window.")
parser.add_argument("--haplotype",action="store_true",help="Calculate D and D+ based on site counts for a haplotype per population.")
parser.add_argument("-haplotype_index", type=str, nargs='?', dest="haplotype_index", help="Haplotype index corresponding to introgression maps.")
parser.add_argument("-p1",type=str,nargs='?',dest="pop_one_name",help="Name of individual in Population 1.")
parser.add_argument("-p2",type=str,nargs='?',dest="pop_two_name",help="Name of individual in Population 2.")
parser.add_argument("-p3",type=str,nargs='?',dest="pop_three_name",help="Name of individual in Population 3.")
parser.add_argument("--d_statistic",action="store_true",help="To run the D statistic add option.")
parser.add_argument("--dplus_statistic",action="store_true",help="To run the D+ statistic add option.")
parser.add_argument("--fd_statistic",action="store_true",help="To run the fD statistic add option.")
parser.add_argument("--fdm_statistic",action="store_true",help="To run the fDM statistic add option.")
parser.add_argument("--df_statistic",action="store_true",help="To run the df statistic add option.")

parser.add_argument("--verbose",action="store_true",help="Verbose mode.")
parser.add_argument("--position",type=str,dest="first_last_position_chromosome_file",help="Where is the file with the first and last position in the vcf for 1000 Genomes Project.")
args=parser.parse_args()

args.chromosome=22
args.window_size=50000

args.names_string="p1[NA06989.22.0]p2[HG00148]p3[Vindija33.19]"

args.first_last_position_chromosome_file="/pl/active/villanea_lab/nickhart_data/dplus_test/first_last_position.csv"
args.human_vcf="/pl/active/villanea_lab/data/modern_data/1000_genomes/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(args.chromosome)
args.archaic_vcf="/pl/active/villanea_lab/data/archaic_data/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr{}_mq25_mapab100.vcf.gz".format(args.chromosome)
args.introgressed_infile="/pl/active/villanea_lab/data/modern_data/steinrucken_neanderthal_SNP_calls/CEU_calls/CEU_lax_chr{}_out.txt"
#outfile="/pl/active/villanea_lab/nickhart_data/dplus_test/p1-${pop_one_name}_p2-${pop_two_name}_p3-${pop_three_name}-2013_window-${window_size}_${chromosome}.csv"

args.d_statistic=True
args.dplus_statistic=True
args.haplotype=True



#Get first and last position of chromosome
first_last_position_chromosome_dict=get_first_last_position_per_chromosome(args.first_last_position_chromosome_file)
first = first_last_position_chromosome_dict[args.chromosome][0]
last  = first_last_position_chromosome_dict[args.chromosome][1]
#Skip the first window and the last window
start = first+args.window_size + 1
stop=start+args.window_size

#What statistics to run
statistics_list=[]
#Do D
if args.d_statistic:
  statistics_list.append("D")
#Do D+
if args.dplus_statistic:
  statistics_list.append("D+")
#Do fD
if args.fd_statistic:
  statistics_list.append("fD")
#Do fDM
if args.fdm_statistic:
  statistics_list.append("fDM")
#Do df
if args.df_statistic:
  statistics_list.append("df")


#Allele sites
allele_sites_list=["ABBA","BABA","BAAA","ABAA"]

#Introgressed sites
if args.haplotype:
  haplotype_string="_"
 haplotype_indexes=[int(haplotype_index) for haplotype_index in args.haplotype_index.split(haplotype_string)[:2]]
  introgressed_tracts={}
  introgressed_tracts[0]=get_introgressed_region_for_haplotype(file_name=args.introgressed_infile,individual=[haplotype_indexes[0]])
  introgressed_tracts[1]=get_introgressed_region_for_haplotype(file_name=args.introgressed_infile,individual=[haplotype_indexes[1]])
else:
  introgressed_tracts=get_introgressed_region_for_frequencies(file_name=args.introgressed_infile,chromosome=args.chromosome)

#Open files with tabix
tb = tabix.open(args.human_vcf)
atb = tabix.open(args.archaic_vcf)

#Open outfile
outfile_delim="\t"
if os.path.exists(args.outfile):
  fout=open(args.outfile,"a+")
else:
  fout=open(args.outfile,"w")
  header=outfile_delim.join(["chromosome","start","stop","number_of_snps","haplotypes_of_p1/p2","introgressed_bases","introgressed_percentage"]+allele_sites_list+statistics_list)
  fout.write("{}\n".format(header))

#Necessary parameters
info_index=7;position_index=1;reference_allele_index=3;alternative_allele_index=4
nucleotides=["A","T","C","G"]
human_individual_index_dict=get_index_populations_dict(vcf_name=args.human_vcf,string=args.names_string,populations=["p1","p2"])
archaic_individual_index_dict=get_index_populations_dict(vcf_name=args.archaic_vcf,string=args.names_string,populations=["p3"])
haplotype_combinations={"second/first":(1,0),"second/second":(1,1)}  #{"first/first":(0,0),"first/second":(0,1)} #,"second/first":(1,0),"second/second":(1,1)}

#Window intervals
#for interval in window_intervals:
while (stop < (last - args.window_size)):
  data=defaultdict(lambda: defaultdict(None))
  allele_sites=defaultdict(lambda: defaultdict(int))
  results=defaultdict(None)
  stats_components=defaultdict(int)
  #Get iterators for human and archaic vcfs
  human_iterator = tb.query(args.chromosome,start,stop)
  archaic_iterator = atb.query(args.chromosome,start,stop)
  for vcf_row in human_iterator:
    #Only SNPs
    if not (("VT=SNP" in vcf_row[info_index])&("AA=" in vcf_row[info_index])):
      continue
    #Get ancestral (P4) and human alleles (P1 and P2)
    ancestral_allele=get_ancestral_allele(vcf_row[position_index],vcf_row[info_index])
    human_reference=vcf_row[reference_allele_index]
    human_alternative=vcf_row[alternative_allele_index]
    #Check biallelic for humans and ancestral allele
    if not ((nucleotides.count(ancestral_allele)==1)&(nucleotides.count(human_reference)==1)&(nucleotides.count(human_alternative)==1)):
      continue
    #Get data for ancestral and humans
    position=int(vcf_row[position_index])
    vcf_row=np.asarray(vcf_row)
    data[position]["ancestral_allele"]=ancestral_allele
    data[position]["human_reference"]=human_reference
    data[position]["human_alternative"]=human_alternative
    if args.haplotype:
      data[position]["pop_one"]=vcf_row[human_individual_index_dict["pop1"]]
      data[position]["pop_two"]=vcf_row[human_individual_index_dict["pop2"]]
    else:
      
data[position]["pop_one_derived_freq"] =get_derived_freq_per_population(array=vcf_row[human_individual_index_dict["pop1"]]
                    ,ancestral_allele=ancestral_allele,reference_allele=human_reference,alternative_allele=human_alternative)
      data[position]["pop_two_derived_freq"]=get_derived_freq_per_population(array=vcf_row[human_individual_index_dict["pop2"]]
                    ,ancestral_allele=ancestral_allele,reference_allele=human_reference,alternative_allele=human_alternative)
  #Positions for all four populations
  positions_list=[]
  #Get data for archaic for altai nean
  for vcf_row in archaic_iterator:
    position=int(vcf_row[position_index])
    #Only SNPs with data from humans and outgroup
    if not position in list(data.keys()):
     continue
    #Reference and alternative alleles for archaic human
    archaic_reference=vcf_row[reference_allele_index]
    archaic_alternative=vcf_row[alternative_allele_index].replace('.',archaic_reference)
    #Is SNP biallelic for all four populations?
    alleles=set([data[position]["ancestral_allele"],data[position]["human_reference"],data[position]["human_alternative"],archaic_reference,archaic_alternative])
    if not len(alleles)==2:
      continue
    positions_list.append(position)
    vcf_row=np.asarray(vcf_row)
    data[position]["archaic_reference"]=archaic_reference
    data[position]["archaic_alternative"]=archaic_alternative
    if args.haplotype:
      data[position]["pop_three"]=vcf_row[archaic_individual_index_dict["pop3"]]
    else:
      data[position]["pop_three_derived_freq"]=get_derived_freq_per_population(array=vcf_row[archaic_individual_index_dict["pop3"]]
                    ,ancestral_allele=data[position]["ancestral_allele"],reference_allele=archaic_reference,alternative_allele=archaic_alternative)
  #Compute statistics for window
  for position in positions_list:
    if args.haplotype:
      #Get site patterns
      for combo in haplotype_combinations.keys():
        pop_one_allele_integer=data[position]["pop_one"].split("|")[haplotype_combinations[combo][0]]
        pop_one_allele=get_allele_from_individual(reference_allele=data[position]["human_reference"]
                                  ,alternative_allele=data[position]["human_alternative"]
                                                    ,allele_integer=pop_one_allele_integer)
        pop_two_allele_integer=data[position]["pop_two"].split("|")[haplotype_combinations[combo][1]]
        pop_two_allele=get_allele_from_individual(reference_allele=data[position]["human_reference"]
                                   ,alternative_allele=data[position]["human_alternative"] 
                                                    ,allele_integer=pop_two_allele_integer)
        pop_three_allele_integer=random.choice(data[position]["pop_three"].split(":")[0].split("/"))
        pop_three_allele=get_allele_from_individual(reference_allele=data[position]["archaic_reference"]
                             ,alternative_allele=data[position]["archaic_alternative"]
                                                    ,allele_integer=pop_three_allele_integer)
        site_pattern = allele_site_pattern(pop1_allele=pop_one_allele,pop2_allele=pop_two_allele
                               ,pop3_allele=pop_three_allele,pop4_allele=data[position]["ancestral_allele"])
        allele_sites[combo][site_pattern]+=1
    else:
      intro_dict=update_introgression_stats_sums(p_one=pop_derived_freq["pop1"],p_two=pop_derived_freq["pop2"]
                                           ,p_three=pop_derived_freq["pop3"],p_four=pop_derived_freq["outgroup"]
                                           ,list=statistics_list)
      stats_components={k: stats_components.get(k,0)+intro_dict.get(k,0) for k in intro_dict.keys()}
  #Next window if window empty
  if len(positions_list)==0:
    start= stop + 1; stop = start+args.window_size
    continue
  #Output results to outfile
  
outline = outfile_delim.join([args.chromosome, str(start), str(stop), str(len(positions_list)), str(introgressed_bases), str(introgressed_bases/args.window_size)])
  if args.haplotype:
    for combo in haplotype_combinations:
      results[combo]=calculate_introgression_stats(allele_sites[combo],statistics_list)
      #Get introgressed bases per window
      introgressed_bases=np.sum([number_intro_bases(window_start,window_stop,intro_tracts[0],intro_tracts[1]) for intro_tracts in introgressed_tracts[haplotype_combinations[combo][1]]])
      outline+="{}{}".format(outfile_delim, outfile_delim.join([combo, str(introgressed_bases), str(introgressed_bases / args.window_size)]))
      outline+="{}{}".format(outfile_delim,outfile_delim.join([str(allele_sites[combo][allele_site]) for allele_site in allele_sites_list]+[str(results[combo][stat]) for stat in statistics_list]))
  else:
    #Calculate the stats for the window
    intro_stats=calculate_introgression_stats(stats_components,statistics_list)
    #Get introgressed bases per window
    introgressed_bases=np.sum([number_intro_bases(window_start,window_stop,intro_tracts[0],intro_tracts[1]) for intro_tracts in introgression_tracts])
    outline += "{}{}".format(outfile_delim, outfile_delim.join([str(introgressed_bases), str(introgressed_bases / args.window_size)]))
    outline+="{}{}".format(outfile_delim,outfile_delim.join([str(intro_stats[key]) for key in statistics_list]))
  fout.write("{}\n".format(outline))
  #Update start and stop positions for new interval
  start= stop + 1
  stop = start + args.window_size
  if args.verbose:
    print(start,stop,len(positions_list))

#Close outfile
fout.close()

