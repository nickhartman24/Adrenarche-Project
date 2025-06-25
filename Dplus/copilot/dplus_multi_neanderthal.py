import argparse
import gzip
import numpy as np
import os
import re
import sys
from collections import defaultdict
import tabix

def print_progress(msg):
    print(f"[INFO] {msg}")

def parse_gene_regions(region_file):
    regions = []
    with open(region_file) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                chrom, start, end = line.strip().split()[:3]
                regions.append((chrom, int(start), int(end)))
    print_progress(f"Loaded {len(regions)} regions from {region_file}")
    return regions

def load_population_samples(metadata_file, pop_codes):
    pop_to_samples = {pop: [] for pop in pop_codes}
    with open(metadata_file) as f:
        header = f.readline().strip().split()
        sample_col = header.index("sample_id")
        pop_col = header.index("pop")
        for line in f:
            cols = line.strip().split()
            if len(cols) <= max(sample_col, pop_col):
                continue
            sample, pop = cols[sample_col], cols[pop_col]
            if pop in pop_to_samples:
                pop_to_samples[pop].append(sample)
    print_progress(f"Sample assignments: {pop_to_samples}")
    return pop_to_samples

def get_population_indices(vcf_file, pop_names):
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split()
                break
    indices = []
    for name in pop_names:
        if name in header:
            indices.append(header.index(name))
        else:
            print_progress(f"WARNING: Sample {name} not found in VCF header {vcf_file}")
    print_progress(f"Population indices for {pop_names}: {indices}")
    return indices

def get_ancestral_allele(info_string):
    try:
        m = re.search(r'AA=([ATCG])', info_string)
        if m:
            return m.group(1)
    except Exception as e:
        print_progress(f"Error getting ancestral allele: {e}")
    return None

def get_derived_allele(ancestral_allele, ref, alt):
    # Returns 0 if ref is derived, 1 if alt is derived
    if ancestral_allele == ref:
        return 1
    elif ancestral_allele == alt:
        return 0
    else:
        return None

def get_derived_allele_freq(genos, ancestral_allele, ref, alt):
    # genos: list of genotype strings (e.g. ['0|1', '1|1', ...])
    derived = get_derived_allele(ancestral_allele, ref, alt)
    if derived is None:
        return None
    count = 0
    total = 0
    for g in genos:
        alleles = re.split('[|/]', g)
        for a in alleles:
            if a == '0':
                if derived == 0:
                    count += 1
            elif a == '1':
                if derived == 1:
                    count += 1
            total += 1
    return count / total if total > 0 else None

def calculate_D_stats(ABBA, BABA, BAAA, ABAA, statistics_list):
    results = {}
    if "D" in statistics_list:
        denom = ABBA + BABA
        results["D"] = (ABBA - BABA) / denom if denom else np.nan
    if "D+" in statistics_list:
        denom = ABBA + BABA + BAAA + ABAA
        results["D+"] = ((ABBA - BABA) + (BAAA - ABAA)) / denom if denom else np.nan
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", type=str, required=True, help="BED file with regions: chrom start end")
    parser.add_argument("--samples_metadata", type=str, required=True, help="Sample metadata")
    parser.add_argument("--pop1_code", type=str, required=True, help="Pop1 code (e.g. GBR)")
    parser.add_argument("--pop2_code", type=str, required=True)
    parser.add_argument("--neanderthal_list", type=str, required=True,
                        help="Comma-separated list of pop3_code:vcf_path, e.g. 'Vindija33.19:/path/to/Vindija/chr{}.vcf.gz,AltaiNea:/path/to/Altai/chr{}.vcf.gz'")
    parser.add_argument("--human_vcf", type=str, required=True, help="VCF path with {} as chrom placeholder")
    parser.add_argument("--outfile", type=str, required=True)
    parser.add_argument("--d_statistic", action="store_true")
    parser.add_argument("--dplus_statistic", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # Parse gene regions
    regions = parse_gene_regions(args.regions)

    # Parse Neanderthal meta-info
    neanderthals = []
    for entry in args.neanderthal_list.split(","):
        code, vcf = entry.split(":")
        neanderthals.append({"code": code, "vcf": vcf})

    # Parse population sample assignments
    # Collect all pop3_codes for Neanderthals
    all_pop3_codes = [n["code"] for n in neanderthals]
    pop_dict = load_population_samples(args.samples_metadata, [args.pop1_code, args.pop2_code] + all_pop3_codes)
    pop1_inds = pop_dict[args.pop1_code]
    pop2_inds = pop_dict[args.pop2_code]
    nean_inds_dict = {n["code"]: pop_dict[n["code"]] for n in neanderthals}

    if len(pop1_inds) == 0 or len(pop2_inds) == 0 or any(len(nean_inds_dict[n["code"]]) == 0 for n in neanderthals):
        print_progress("ERROR: One or more populations have no individuals specified.")
        sys.exit(1)
    # Only use the first sample for each Neanderthal
    for n in neanderthals:
        if len(nean_inds_dict[n["code"]]) > 1:
            print_progress(f"WARNING: Multiple samples found for {n['code']}; using only the first.")
            nean_inds_dict[n["code"]] = nean_inds_dict[n["code"]][:1]

    # Setup statistics
    statistics_list = []
    if args.d_statistic:
        statistics_list.append("D")
    if args.dplus_statistic:
        statistics_list.append("D+")

    # Open output
    fout = open(args.outfile, "w")
    # Output header: region info + n_snps + D/D+ for each Neanderthal
    header = ["chrom","start","end","n_snps"]
    for n in neanderthals:
        for stat in statistics_list:
            header.append(f"{stat}_{n['code']}")
    fout.write("\t".join(header) + "\n")

    # Loop over regions
    for region_number, (chrom, start, end) in enumerate(regions):
        print_progress(f"Processing region {region_number+1}/{len(regions)}: {chrom}:{start}-{end}")

        # Format per-chromosome human VCF file path
        human_vcf = args.human_vcf.format(chrom)
        try:
            tb = tabix.open(human_vcf)
        except Exception as e:
            print_progress(f"Tabix open failed for human VCF {human_vcf}: {e}")
            continue

        # Get sample indices for pop1/pop2 for this file
        pop1_idx = get_population_indices(human_vcf, pop1_inds)
        pop2_idx = get_population_indices(human_vcf, pop2_inds)

        if len(pop1_idx) == 0 or len(pop2_idx) == 0:
            print_progress("ERROR: pop1 or pop2 has no valid sample indices in the VCF header.")
            continue

        # Query region in human VCF
        try:
            human_iter = tb.query(str(chrom), start, end)
        except Exception as e:
            print_progress(f"Tabix query failed for {chrom}:{start}-{end}: {e}")
            continue

        # SNPs per region
        region_data = {}
        for vcf_row in human_iter:
            if not vcf_row[0] == str(chrom):
                continue
            info = vcf_row[7]
            if not ("VT=SNP" in info and "AA=" in info):
                continue
            ancestral_allele = get_ancestral_allele(info)
            ref = vcf_row[3]
            alt = vcf_row[4]
            if not (ancestral_allele and ref in "ATCG" and alt in "ATCG"):
                continue
            pos = int(vcf_row[1])
            # Extract genotypes for each pop
            geno_pop1 = [vcf_row[i] for i in pop1_idx]
            geno_pop2 = [vcf_row[i] for i in pop2_idx]
            region_data[pos] = {
                "ancestral": ancestral_allele,
                "ref": ref,
                "alt": alt,
                "geno_pop1": geno_pop1,
                "geno_pop2": geno_pop2
            }
        print_progress(f"  Loaded {len(region_data)} human SNPs for region.")

        snp_count = len(region_data)
        nean_stats = {}

        # For each Neanderthal, get indices, load archaic data, and compute D-stats
        for n in neanderthals:
            nean_vcf = n["vcf"].format(chrom)
            try:
                atb = tabix.open(nean_vcf)
            except Exception as e:
                print_progress(f"Tabix open failed for {n['code']} VCF {nean_vcf}: {e}")
                nean_stats[n["code"]] = {stat: "nan" for stat in statistics_list}
                continue

            pop3_inds = nean_inds_dict[n["code"]]
            pop3_idx = get_population_indices(nean_vcf, pop3_inds)
            if len(pop3_idx) == 0:
                print_progress(f"ERROR: {n['code']} has no valid sample indices in VCF header.")
                nean_stats[n["code"]] = {stat: "nan" for stat in statistics_list}
                continue

            try:
                archaic_iter = atb.query(str(chrom), start, end)
            except Exception as e:
                print_progress(f"Tabix query failed for {n['code']} {chrom}:{start}-{end}: {e}")
                nean_stats[n["code"]] = {stat: "nan" for stat in statistics_list}
                continue

            # Load archaic data and match to human positions
            archaic_pos_set = set()
            for vcf_row in archaic_iter:
                if not vcf_row[0] == str(chrom):
                    continue
                pos = int(vcf_row[1])
                if pos not in region_data:
                    continue
                ref = vcf_row[3]
                alt = vcf_row[4].replace('.', ref)
                geno_pop3 = [vcf_row[i] for i in pop3_idx]
                region_data[pos][f"ref3_{n['code']}"] = ref
                region_data[pos][f"alt3_{n['code']}"] = alt
                region_data[pos][f"geno_pop3_{n['code']}"] = geno_pop3
                archaic_pos_set.add(pos)
            print_progress(f"  {n['code']}: Loaded archaic data for {len(archaic_pos_set)} SNPs.")

            # Compute stats
            ABBA = BABA = BAAA = ABAA = 0
            for pos in archaic_pos_set:
                d = region_data[pos]
                try:
                    anc = d["ancestral"]
                    ref = d["ref"]
                    alt = d["alt"]
                    ref3 = d[f"ref3_{n['code']}"]
                    alt3 = d[f"alt3_{n['code']}"]
                    # Calculate derived allele frequency for each pop at this SNP
                    p1_freq = get_derived_allele_freq(d["geno_pop1"], anc, ref, alt)
                    p2_freq = get_derived_allele_freq(d["geno_pop2"], anc, ref, alt)
                    p3_freq = get_derived_allele_freq(d[f"geno_pop3_{n['code']}"], anc, ref3, alt3)
                    if None in (p1_freq, p2_freq, p3_freq):
                        continue
                    # Use 0.5 as threshold for presence/absence
                    p1 = int(p1_freq > 0.5)
                    p2 = int(p2_freq > 0.5)
                    p3 = int(p3_freq > 0.5)
                    # Patterns
                    if p3 == p2 and p1 != p2:
                        ABBA += 1
                    elif p3 == p1 and p1 != p2:
                        BABA += 1
                    # D+ extra patterns (optional, simple version)
                    if p3 == p2 and p1 != p2:
                        BAAA += 1
                    elif p3 == p1 and p1 != p2:
                        ABAA += 1
                except Exception as e:
                    print_progress(f"Error calculating stats for {n['code']} at pos {pos}: {e}")

            stats = calculate_D_stats(ABBA, BABA, BAAA, ABAA, statistics_list)
            nean_stats[n["code"]] = stats

        # Output
        row = [chrom, start, end, snp_count]
        for n in neanderthals:
            for stat in statistics_list:
                val = nean_stats[n["code"]].get(stat, "nan")
                row.append(val)
        fout.write("\t".join(map(str, row)) + "\n")

    fout.close()
    print_progress(f"All done! Output written to {args.outfile}")

if __name__ == "__main__":
    main()
