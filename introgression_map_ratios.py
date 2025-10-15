import csv

# User-defined variables
gene = "CYB5A"
pop = "CEU"
chromosome = "chr18"
gene_first_position = 71920750
gene_last_position = 71959250

input_file = f"/pl/active/villanea_lab/data/modern_data/steinrucken_neanderthal_SNP_calls/{pop}_calls/{pop}_lax_{chromosome}_out.txt"
output_file = f"{chromosome}_{pop}_ratio_output.txt"
filtered_file = f"{gene}_{pop}_ratio_output_filtered.txt"

# Step 1: Convert individual probabilities to 1 (likely introgressed) and 0 (unlikely introgressed), with threshold of 0.89. 
#Then create an output file that contains the percentage of “1” individuals out of the whole pop sample. 
with open(input_file, newline='') as infile, open(output_file, "w", newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    header = next(reader)
    n = len(header) - 2  # Number of individuals

    # Write new header
    writer.writerow([header[0], header[1], "ratio"])

    for row in reader:
        chrom = row[0]
        pos = row[1]
        ones = 0
        for val in row[2:]:
            try:
                num = float(val)
                if num >= 0.89:
                    ones += 1
            except ValueError:
                continue  # skip non-numeric
        ratio = ones / n if n > 0 else 0
        writer.writerow([chrom, pos, "%.6f" % ratio])

# Step 2: Filter the ratio output of the whole chromosome down to just the gene positions. This allows for quick reference to the gene, while still having the previous file for comparison to whole chromosome. 
with open(output_file, newline='') as infile, open(filtered_file, "w", newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    header = next(reader)
    writer.writerow(header)

    for row in reader:
        pos = int(row[1])
        if gene_first_position <= pos <= gene_last_position:
            writer.writerow(row)

print(f"Done! Outputs written to:\n- {output_file}\n- {filtered_file}")

