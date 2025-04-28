# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:22:10 2025

@author: George_K
"""

import pysam
import pandas as pd
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process a BAM file and split reads by region based on spot barcodes.")
parser.add_argument("--data_wd", default=None, required=False, help="Path to the directory containing the BAM and CSV files.")
parser.add_argument("--bam_file", required=True, help="Name of the BAM file")
parser.add_argument("--region_file", required=True, help="Name of the CSV file - region annotation")
parser.add_argument("--output_dir",required=True, help="Path to the output directory of the new BAM files.")

args = parser.parse_args()


# Path to the data
data_wd = args.data_wd
output_dir = args.output_dir

if (data_wd is not None):
    bam_file = data_wd + args.bam_file  # Full path to the BAM file
    region_file = data_wd + args.region_file  # Full path to the CSV file
else:
    bam_file = args.bam_file
    region_file = args.region_file

# Path to the data
# data_wd = '/home/georgios/NeoSpan/Data/'

# Load the BAM file
# bam_file = data_wd + "Visium_FFPE_Human_Prostate_Cancer_possorted_genome_bam.bam"  # Path to your BAM file
samfile = pysam.AlignmentFile(bam_file, "rb")


# Load the CSV file containing the barcode to region mapping
# csv_file = data_wd + "Pathology.csv"  # Path to your CSV file
barcode_region_df = pd.read_csv(region_file, delimiter = ';')

# Get unique regions from the DataFrame
regions = barcode_region_df['Pathology'].unique()


# Create a dictionary to map each region to its output file
region_to_file = {}

# Open a file for each region
for region in regions:
    output_filename = output_dir + f"{region}_output.bam"
    region_to_file[region] = pysam.AlignmentFile(output_filename, "wb", template=samfile)

# Loop through the reads in the BAM file
for read in samfile.fetch():
    if read.has_tag("CB"):
        barcode = read.get_tag("CB")  # Extract barcode from the read's CB tag
        
        # Find the region for this barcode (assuming barcode_region_df has 'Barcode' and 'Pathology' columns)
        region = barcode_region_df.loc[barcode_region_df['Barcode'] == barcode, 'Pathology']
        
        # Write the read to the corresponding region's file
        if not region.empty:
            region = region.values[0]
            region_to_file[region].write(read)
        
            

# Close all output files
for region_f in region_to_file.values():
    region_f.close()


# # Prepare to write out two separate BAM files
# normal_bam = data_wd + "normal_reads.bam"
# tumor_bam = data_wd + "tumor_reads.bam"

# # Create new BAM files (as temporary files)
# normal_outfile = pysam.AlignmentFile(normal_bam, "wb", header=samfile.header)
# tumor_outfile = pysam.AlignmentFile(tumor_bam, "wb", header=samfile.header)

# # Create a set of barcodes that correspond to normal and tumor regions
# normal_barcodes = set(barcode_region_df[barcode_region_df['Pathology'] == 'Normal gland']['Barcode'])
# tumor_barcodes = set(barcode_region_df[barcode_region_df['Pathology'] != 'Normal gland']['Barcode'])

# # Loop through the reads in the BAM file
# for read in samfile.fetch():
    
#     if read.has_tag("CB"):
#         barcode = read.get_tag("CB")  # Extract barcode from the read's CB tag
#         if barcode in normal_barcodes:
#             normal_outfile.write(read)
#         elif barcode in tumor_barcodes:
#             tumor_outfile.write(read)

# Close all files
samfile.close()
# normal_outfile.close()
# tumor_outfile.close()

print("BAM files for the regions have been created.")

