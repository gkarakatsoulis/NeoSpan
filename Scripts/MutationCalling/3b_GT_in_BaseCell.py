#!/usr/bin/env python3

import argparse
import os

def parse_base_counts(bc):
    """
    Parse the BC field into counts for A, C, T, G, I, D, N, O.
    If the BC field has fewer than 8 values, pad it with zeros.
    """
    counts = list(map(int, bc.split(':')))
    while len(counts) < 8:
        counts.append(0)
    return counts

def calculate_gt(ref, bc):
    """
    Calculate the genotype (GT) based on reference allele and base counts.
    """
    a_count, c_count, t_count, g_count, i_count, d_count, n_count, o_count = parse_base_counts(bc)
    
    base_counts = {'A': a_count, 'C': c_count, 'T': t_count, 'G': g_count}
    dp = sum(base_counts.values())  # Total depth
    
    if dp == 0:
        return "./.", "."

    if ref in base_counts:
        del base_counts[ref]
    
    alt_allele = max(base_counts, key=base_counts.get)
    alt_count = base_counts[alt_allele]
    
    alt_ratio = alt_count / dp if dp > 0 else 0
    
    if alt_ratio > 0.8:
        gt = '1/1'
    elif alt_ratio > 0:
        gt = '0/1'
    else:
        gt = '0/0'
        alt_allele = '.'



    # if alt_ratio > 0.8:
    #     gt = '1/1'
    # elif alt_ratio > 0.2:
    #     gt = '0/1'
    # else:
    #     gt = '0/0'
    #     alt_allele = '.'
    
    return gt, alt_allele

def process_base_count_matrix(input_file, output_file):
    """
    Process the base count matrix, calculate GT, and write to output TSV.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Adjusted header (removed DP, moved GT inside INFO)
        outfile.write("#CHROM\tPOS\tREF\tALT\tINFO\tCounts\n")

        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, info, counts = fields[0], fields[1], fields[2], fields[3], fields[4]
            
            counts_fields = counts.split('|')
            if len(counts_fields) < 4:
                print(f"Skipping line due to missing INFO fields: {line.strip()}")
                continue
            bc = counts_fields[3]

            try:
                gt, alt_allele = calculate_gt(ref, bc)
            except ValueError:
                print(f"Skipping line due to invalid BC field: {line.strip()}")
                continue
            
            # Modify INFO column by prepending GT
            updated_info = f"GT|{info}"
            updated_counts = f"{gt}|{counts}"

            outfile.write(f"{chrom}\t{pos}\t{ref}\t{alt_allele}\t{updated_info}\t{updated_counts}\n")

def process_folder(input_folder):
    """
    Process all files in the input folder, adding GT field and saving new files with '_GT' suffix.
    """
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        return
    
    files = [f for f in os.listdir(input_folder) if f.endswith(".txt") or f.endswith(".tsv")]
    
    if not files:
        print(f"No .txt or .tsv files found in {input_folder}.")
        return
    
    output_folder = input_folder + 'GT'
    os.makedirs(output_folder, exist_ok=True)

    for filename in files:
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename.rsplit('.', 1)[0] + "_GT.tsv")

        print(f"Processing {input_path} -> {output_path}")
        process_base_count_matrix(input_path, output_path)
    
    print("All files processed.")

def main():
    parser = argparse.ArgumentParser(description="Modify base count matrix files: remove DP, move GT under INFO.")
    parser.add_argument("--indir", required=True, help="Path to the input folder containing base count matrix files.")
    args = parser.parse_args()
    
    process_folder(args.indir)

if __name__ == "__main__":
    main()
