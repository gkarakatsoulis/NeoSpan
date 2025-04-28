import pandas as pd
import argparse

def feat_spot_create(infile, file_dir=None, sep='\t', outfile="feature_spot_matrix.tsv"):

    # Construct the full path to the input file
    if file_dir is not None:
        infile = file_dir + infile
    
    # Identify the first line with actual column headers
    with open(infile, "r") as file:
        for line_number, line in enumerate(file):
            if line.startswith("#CHROM"):  # This is the actual header line
                header_line = line_number
                break
            
    
    
    # Read the input file
    df = pd.read_csv(infile, sep=sep, skiprows=header_line)
    

    
    print(df.columns)
    
    # Create a unique variant identifier
    df["Variant_ID"] = df["#CHROM"] + "_" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]

    # Pivot to create the feature-spot matrix
    feature_spot_matrix = df.pivot(index="Variant_ID", columns="Cell_types", values="VAF")

    # Fill missing values with 0 (if a spot has no detected mutation)
    feature_spot_matrix = feature_spot_matrix.fillna(0)

    # Save the feature-spot matrix
    feature_spot_matrix.to_csv(outfile, sep="\t")

    # # Display the final matrix
    # print("Feature-Spot Matrix:")
    # print(feature_spot_matrix.head())

def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to construct the feature-spot matrix')
    parser.add_argument('--infile', type=str, help='Input file with the mutations and spot information', required=True)
    parser.add_argument('--file_dir', type=str, help='Path for the input file', required=False)
    parser.add_argument('--outfile', type=str, help='Output file prefix', required=False)
    return parser

def main():
    # Initialize the argument parser
    parser = initialize_parser()
    args = parser.parse_args()

    # Call the function to create the feature-spot matrix
    feat_spot_create(
        infile=args.infile,
        file_dir=args.file_dir,
        outfile=args.outfile
    )

if __name__ == "__main__":
    main()
