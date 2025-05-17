
#!/usr/bin/env python3
import sys
import argparse
import csv

VCF_HEADER = """##fileformat=VCFv4.0
##source=VarscanSomatic
##reference=ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz
##phasing=none
##center=genome.wustl.edu
##FILTER=<ID=PASS,Description="Passed all filters">
##FILTER=<ID=VarscanHighConfidenceIndel,Description="Filter description">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles 0/1/2/3...">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|...">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb
"""

def parse_args():
    parser = argparse.ArgumentParser(description="Convert custom TSV to VCF")
    parser.add_argument("input_file", help="Input TSV file")
    parser.add_argument("output_file", help="Output VCF file")
    return parser.parse_args()

def process_line(row):
    chrom = row['#CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    if alt == '.':
        alt = row['REF']

    filt = row['FILTER'].split(',')[0] if row['FILTER'] else 'PASS'
    bam_gt = row['b']
    info_dict = dict(zip(row['INFO'].split('|'), bam_gt.split('|')))
    dp = info_dict.get('DP', '.')
    bc = info_dict.get('BC', '0:0:0:0:0:0:0:0')
    ad_list = bc.split(':')[:4]  # A:C:T:G
    ad = ','.join(ad_list)
    gt = info_dict.get('GT', './.')

    return {
        'CHROM': chrom,
        'POS': pos,
        'ID': '.',
        'REF': ref,
        'ALT': alt,
        'QUAL': '.',
        'FILTER': filt,
        'INFO': '.',
        'FORMAT': 'GT:DP:AD',
        'SAMPLE': f"{gt}:{dp}:{ad}"
    }

def main():
    args = parse_args()

    with open(args.input_file) as tsvfile, open(args.output_file, 'w') as vcffile:
        lines = tsvfile.readlines()

        # Skip meta-information and extract header
        header_line = next((line for line in lines if line.startswith("#CHROM")), None)
        if not header_line:
            raise ValueError("No header line starting with '#CHROM' found in the TSV file.")

        headers = header_line.strip().split("\t")
        reader = csv.DictReader(lines[lines.index(header_line)+1:], fieldnames=headers, delimiter="\t")

        vcffile.write(VCF_HEADER)

        for row in reader:
            if row.get('ALT', '.') == '.':
                row['ALT'] = row['REF']

            record = process_line(row)
            vcf_line = "\t".join([
                record['CHROM'],
                record['POS'],
                record['ID'],
                record['REF'],
                record['ALT'],
                record['QUAL'],
                record['FILTER'],
                record['INFO'],
                record['FORMAT'],
                record['SAMPLE']
            ])
            vcffile.write(vcf_line + "\n")


if __name__ == "__main__":
    main()





# import pandas as pd
# import argparse

# def convert_scomatic_tsv_to_vcf(infile: str, outfile: str):
#     """
#     Convert a SComatic TSV-like output file to a VCF format with proper genotype columns.
#     Modified version for pVACseq compatibility.
#     """
#     print("Reading input file...")
    
#     # Read metadata lines (lines starting with ##)
#     metadata_lines = []
#     with open(infile, "r") as file:
#         for line in file:
#             if line.startswith("##"):
#                 metadata_lines.append(line.strip())
#             elif line.startswith("#CHROM"):  # This is the header line
#                 header_line = line.strip().split("\t")
#                 break
    
#     # Read the actual TSV data (skip metadata)
#     df = pd.read_csv(infile, sep="\t", comment="#", names=header_line)
    
#     print("Processing data...")
    
#     # Define VCF header with proper FORMAT definitions
#     vcf_header = "##fileformat=VCFv4.2\n"
#     vcf_header += "##source=SComatic\n"
#     vcf_header += "\n".join(metadata_lines) + "\n"
    
#     # Add required FORMAT definitions for pVACseq
#     vcf_header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
#     vcf_header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
#     vcf_header += "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
    
#     # Get sample names (assuming they're all columns after INFO)
#     sample_columns = [col for col in df.columns if col not in 
#                      ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
    
#     # Build the header line
#     vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
#     if sample_columns:
#         vcf_header += "\t" + "\t".join(sample_columns)
#     vcf_header += "\n"
    
#     # Convert TSV to VCF format with proper genotype columns
#     vcf_rows = []
#     for _, row in df.iterrows():
#         chrom = row["#CHROM"]
#         pos = row["POS"]
#         ref = row["REF"]
#         alt = row["ALT"]
#         filter_status = row["FILTER"]
        
#         # Construct INFO field (from non-sample, non-standard columns)
#         info_fields = []
#         for col in df.columns:
#             if col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + sample_columns:
#                 continue
#             value = row[col]
#             if pd.notna(value):
#                 info_fields.append(f"{col}={value}")
        
#         info = ";".join(info_fields)
        
#         # Process sample genotypes
#         format_fields = "GT:DP:AD"  # Minimal required fields for pVACseq
#         sample_data = []
        
#         for sample in sample_columns:
#             # Extract genotype information (assuming format: GT|DP|...)
#             if pd.notna(row[sample]):
#                 sample_parts = str(row[sample]).split("|")
#                 gt = sample_parts[0] if len(sample_parts) > 0 else "./."
#                 dp = sample_parts[1] if len(sample_parts) > 1 else "."
#                 ad = "0,0"  # Placeholder - adjust based on your data
                
#                 # For pVACseq, we need at least GT:DP:AD format
#                 sample_data.append(f"{gt}:{dp}:{ad}")
#             else:
#                 sample_data.append("./.:.:.")
        
#         # Build the VCF line
#         vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\t{format_fields}"
#         if sample_data:
#             vcf_line += "\t" + "\t".join(sample_data)
            
#         vcf_rows.append(vcf_line)
    
#     print("Writing VCF file...")
    
#     # Write to VCF file
#     with open(outfile, "w") as vcf:
#         vcf.write(vcf_header)
#         vcf.write("\n".join(vcf_rows) + "\n")
    
#     print(f"VCF file saved as {outfile}")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Convert SComatic TSV to pVACseq-compatible VCF.")
#     parser.add_argument("infile", help="Path to the input TSV file.")
#     parser.add_argument("outfile", help="Path to the output VCF file.")
#     args = parser.parse_args()
    
#     convert_scomatic_tsv_to_vcf(args.infile, args.outfile)


# # import pandas as pd
# # import argparse

# # def convert_scomatic_tsv_to_vcf(infile: str, outfile: str):
# #     """
# #     Convert a SComatic TSV-like output file to a VCF format.

# #     :param infile: Path to the input file.
# #     :param outfile: Path to the output VCF file.
# #     """
# #     print("Reading input file...")
    
# #     # Read metadata lines (lines starting with ##)
# #     metadata_lines = []
# #     with open(infile, "r") as file:
# #         for line in file:
# #             if line.startswith("##"):
# #                 metadata_lines.append(line.strip())
# #             elif line.startswith("#CHROM"):  # This is the header line
# #                 header_line = line.strip().split("\t")
# #                 break
    
# #     # Read the actual TSV data (skip metadata)
# #     df = pd.read_csv(infile, sep="\t", comment="#", names=header_line)
    
# #     # Debug: Print column names to verify
# #     # print("Columns in the TSV file:", df.columns.tolist())
    
# #     print("Processing data...")
    
# #     # Define VCF header
# #     vcf_header = "##fileformat=VCFv4.2\n"
# #     vcf_header += "##source=SComatic\n"
# #     vcf_header += "\n".join(metadata_lines) + "\n"  # Add metadata lines
# #     vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    
# #     # Convert TSV to VCF format
# #     vcf_rows = []
# #     for _, row in df.iterrows():
# #         chrom = row["#CHROM"]
# #         pos = row["POS"]
# #         ref = row["REF"]
# #         alt = row["ALT"]
# #         filter_status = row["FILTER"]
        
# #         # Construct INFO field
# #         info_fields = []
# #         for col in df.columns[6:]:  # Skip the first 6 columns (#CHROM, POS, REF, ALT, ALT_Mut, FILTER)
# #             if col == "INFO":
# #                 continue  # Skip the INFO column if it exists
# #             value = row[col]
# #             if pd.notna(value):  # Only add non-NaN values
# #                 info_fields.append(f"{col}={value}")
        
# #         info = ";".join(info_fields)
        
# #         # Append the VCF row
# #         vcf_rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}")
    
# #     print("Writing VCF file...")
    
# #     # Write to VCF file
# #     with open(outfile, "w") as vcf:
# #         vcf.write(vcf_header)
# #         vcf.write("\n".join(vcf_rows) + "\n")
    
# #     print(f"VCF file saved as {outfile}")

# # if __name__ == "__main__":
# #     parser = argparse.ArgumentParser(description="Convert a SComatic TSV output file to a VCF format.")
# #     parser.add_argument("infile", help="Path to the input TSV file.")
# #     parser.add_argument("outfile", help="Path to the output VCF file.")
# #     args = parser.parse_args()
    
# #     convert_scomatic_tsv_to_vcf(args.infile, args.outfile)
