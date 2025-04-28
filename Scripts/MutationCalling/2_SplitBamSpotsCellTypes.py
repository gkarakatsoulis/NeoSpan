import pysam
import pandas as pd
import argparse
import timeit
import sys
import numpy as np
import os
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


# ============================================================================
# Function: clean_and_append_tissue
# ============================================================================

# Function to clean a column and append tissue information.

def clean_and_append_tissue(df, column, tissue):
    
    df[f'{column}_clean'] = df[column].str.replace(' ', '_', regex=True)
    if tissue is not None:
        tissue = tissue.replace(" ", "_")
        df[f'{column}_clean'] = tissue + '__' + df[f'{column}_clean'].astype(str)
    return df

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


# ============================================================================
# Function: meta_to_dict
# ============================================================================

# Arguments:
# spot_txt: the path for the spot barcodes
# cell_type_txt: the path for the cell type barcodes
# tissue: defines the tissue - can be None-

#SOS: If both spot and cell barcodes are provided, then we also need
## to map each cell to a specific spot.
## This means that we also need a csv file for the mapping above
## Therefore, either all spot_txt, cell_type_txt and spot_cell_map will be non-None
## or spot_txt and cell_type_txt will be mutually exclusive

def meta_to_dict(tissue, spot_txt = None, cell_type_txt = None, spot_cell_map = None):
    
    if (spot_txt is not None) and (cell_type_txt is not None) and (spot_cell_map is None):
        print ('Error: No mapping between the spot and cell barcodes provided')
        sys.exit()
    
    # Initialize DataFrames
    metadata_map, metadata_spot, metadata_cell_type = None, None, None
    
    if spot_cell_map is not None:
         # Read the meta data as a Pandas DataFrame  
         metadata_map = pd.read_csv(spot_cell_map, delimiter = "\t")
         metadata_map = clean_and_append_tissue(metadata_map, 'Spot', tissue)
         metadata_map = clean_and_append_tissue(metadata_map, 'Cell_type', tissue)
         
        
    
    if spot_txt is not None:
         # Read the meta data as a Pandas DataFrame  
         metadata_spot = pd.read_csv(spot_txt, delimiter = ",")
         
         # Clean index column
         metadata_spot['Index_clean_spot'] = metadata_spot['Index'].str.replace('-.*$','',regex=True)
         
         # Clean the Spot column --> Name it as Spot_clean
         metadata_spot = clean_and_append_tissue(metadata_spot, 'Spot', tissue)


    if cell_type_txt is not None:
        # Read the meta data as a Pandas DataFrame  
        metadata_cell_type = pd.read_csv(cell_type_txt, delimiter = "\t")
        
        # Clean index column
        metadata_cell_type['Index_clean_cell_type'] = metadata_cell_type['Index'].str.replace('-.*$','',regex=True)
        
        # Clean the Spot column --> Name it as Spot_clean
        metadata_cell_type = clean_and_append_tissue(metadata_cell_type, 'Cell_type', tissue)
      
    
    if (spot_txt is not None) and (cell_type_txt is not None):
        metadata = pd.merge(metadata_spot, metadata_map, left_on='Spot_clean', right_on='Spot_clean', how='inner')
        metadata = pd.merge(metadata, metadata_cell_type, left_on='Cell_type_clean', right_on='Cell_type_clean', how='inner')
        metadata['Index_clean'] = "Spot_" + metadata['Index_clean_spot'].astype(str) + "_CellType_" + metadata['Index_clean_cell_type'].astype(str)
        metadata['Annotation'] = "Spot_" + metadata['Spot_clean'].astype(str) + "_CellType_" + metadata['Cell_type_clean'].astype(str)
    elif spot_txt is not None:
        metadata = metadata_spot.rename(columns = {'Index_clean_spot' : 'Index_clean', 'Spot_clean' : 'Annotation'})
    elif cell_type_txt is not None:
        metadata = metadata_cell_type.rename(columns = {'Index_clean_cell_type' : 'Index_clean', 'Cell_type_clean' : 'Annotation'})
    
    
    

	# Create dictionary with cell types and cell barcodes
    DICT = metadata.set_index('Index_clean')['Annotation'].to_dict()
    ALL_Annotations = metadata['Annotation'].unique()
    
    del metadata
    
    # Return the dictionary and the unique cell types
    return(DICT, ALL_Annotations) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# ============================================================================
# Function: split_bam
# ============================================================================


# Arguments
# bam_path: Path to the input BAM file.
# outdir: Directory where output BAM files and the report will be saved.
# donor: Identifier for the donor (used in naming output files).
# tissue: Tissue type (used in naming output files and processing metadata).
# max_NM: Maximum number of mismatches allowed in a read.
# max_NH: Maximum number of hits (alignments) allowed for a read.
# min_MAPQ: Minimum mapping quality score required for a read.
# n_trim: Number of bases to trim from the start and end of each read.
# spot_txt: the path for the spot barcodes
# cell_type_txt: the path for the cell type barcodes
# spot_cell_map: the path for mapping spot barcodes to cell barcodes


def split_bam(bam, outdir, donor, tissue, max_NM, max_NH, min_MAPQ, n_trim, spot_txt = None, cell_type_txt = None, spot_cell_map = None):
    
    start = timeit.default_timer()
    
    # 1. Check if BAM split needs to be done (spot and/or cell type)
    # If not, print a message and exit
    if (spot_txt is None) and (cell_type_txt is None):
        print('Message: No spot/cell type files provided for splitting the BAM file.')
        sys.exit()
    
	# 1. Transform table to dictionary
    DICT, ALL_Annotations = meta_to_dict(tissue, spot_txt, cell_type_txt, spot_cell_map)
    
    if len(DICT) < 1:
        print('Warning: No barcodes found in the --meta file.')
        sys.exit()
    
    # 2. Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    
	# 3. Open infile
    infile=pysam.Samfile(bam, "rb") 
    
    if not pysam.AlignmentFile(bam).has_index():
        print(f"Index file not found for {bam}. Creating index...")
        pysam.index(bam)
        print(f"Index file created: {bam}.bai")
        
        # Reopen the BAM file to ensure the index is recognized
        infile.close()  # Close the current file handle
        infile = pysam.Samfile(bam, "rb")  # Reopen the BAM file


	# 4. Create and open out bam files
    DICT_files = {}
    for spot_cell_type in ALL_Annotations:
	    outfile_bam = f"{outdir}/{donor}.{spot_cell_type}.bam"
	    DICT_files[spot_cell_type] = pysam.AlignmentFile(outfile_bam, "wb",template=infile)

	# 5. Start read counts
    total_reads = 0
    FILTER_dict = {'Total_reads' : 0,'Pass_reads' : 0, 'CB_not_found' : 0, 'CB_not_matched' : 0} # To store filter reasons
	
	# 6. Check reads and split them in bam files    
    for read in infile.fetch():
	    FILTER_dict['Total_reads'] += 1 # To count the total number of reads analysed
	    total_reads += 1

	    if (total_reads % 5000000 == 0):
		    print ('Number of reads already processed: ' + str(total_reads))

		# Check if CB tag is present
	    try:
		    barcode = read.opt("CB")
	    except:
	        FILTER_dict["CB_not_found"] += 1
	        continue 


		# Check if CB code matches with an annotated spot and/or cell type
	    barcode = barcode.split("-")[0]
	    try:
	    	SPOT_CELL_TYPE = DICT[barcode]
	    except:
	    	FILTER_dict["CB_not_matched"] += 1
	    	continue 

		# Final filters
	    FILTER = [] 
		# Check number of mismatches in the read
	    if (max_NM is not None):
	    	try:
	    		if (read.opt("nM") > max_NM):
	    			FILTER.append('nM')
	    	except KeyError:
	    		FILTER.append('nM_not_found')
	    # Check number of hits
	    if (max_NH is not None):
	    	try:
	    		if (read.opt("NH") > max_NH):
	    			FILTER.append('NH')
	    	except KeyError:
	    		FILTER.append('NH_not_found')

		# Check mapping quality
	    if (min_MAPQ > 0):
	    	try:
	    		if (read.mapq < min_MAPQ):
	    			FILTER.append('MAPQ')
	    	except:
	    		FILTER.append('MAPQ_not_found')

		# Making a decision about the read (filtered or not)
	    if (len(FILTER) > 0): # If there are reasons to filter read
	    	FILTER_dict[';'.join(FILTER)] = FILTER_dict.get(';'.join(FILTER), 0) + 1
	    	continue
	    else:
	    	FILTER_dict['Pass_reads'] += 1

		# Only for PASS reads
		# Trim last and first bases of the read (reduce quality) if specified
		# It does not consider the soft-clip bases at the beginning/end of the read for the counting
		# It also considers the expected adapter lengths (up to 30) of 10x library prep to remove bases in long soft-clip sequences
	    if (n_trim > 0):
	        CIGAR_tuple = read.cigartuples
	        if (CIGAR_tuple != None and len(CIGAR_tuple) > 1):
	    		# Check the number of bases to trim at the beginning of the read
				# Extension of the soft-clipped bases
	    	    if (CIGAR_tuple[0][0] == 4):
					# As there are some library preparation adapters inside the 10x prep protocol that 
					# are not properly removed (up to 30 bps), if we observed that the softclip bases are around more than 20, we 
					# prefer to be conservative and remove the expected 30 pbs + the n_trim parameter
			        if (CIGAR_tuple[0][1] >= 20 and CIGAR_tuple[0][1] < 30):
			            trim_start = 30 + n_trim
			        else:
			    	    trim_start = CIGAR_tuple[0][1] + n_trim
	    	    else:
	    	        trim_start = n_trim

				# Check the number of bases to trim at the end of the read
				# Extension of the soft-clipped bases
	    	    if (CIGAR_tuple[-1][0] == 4):	
					# As there are some library preparation adapters inside the 10x prep protocol that 
					# are not properly removed (up to 30 bps), if we observed that the softclip bases are around more than 20, we 
					# prefer to be conservative and remove the expected 30 pbs + the n_trim parameter
	    	    	if (CIGAR_tuple[-1][1] >= 20 and CIGAR_tuple[-1][1] < 30):
	    	    		trim_end = 30 + n_trim
	    	    	else:
	    	    		trim_end = CIGAR_tuple[-1][1] + n_trim
	    	    else:
	    	    	trim_end = n_trim
	        else:
	        	trim_start = n_trim
	        	trim_end = n_trim

			# Set base qualities to 0 at the beginning and the end of the read if specified
			# Get first and last qualities
	        Q_indexes = list(np.arange(trim_start)) + list((np.arange(trim_end) + 1)*-1)
	        Q = read.query_qualities
	        for Qi in Q_indexes:
	        	Q[Qi] = 0

			# Substitute qualities
	        read.query_qualities = Q

		# Print passed read
	    DICT_files[SPOT_CELL_TYPE].write(read)

	# 7. Close opened files
    infile.close()
    for key in DICT_files.keys():
    	DICT_files[key].close()

	# 8. Get and create report
    outfile_report = "{}/{}.{}.report.txt".format(outdir, donor, spot_cell_type)
    stop = timeit.default_timer()
    endtime = round((stop - start),2)
    FILTER_dict['Total_time'] = endtime

    data_df = pd.DataFrame([FILTER_dict])
    data_df.to_csv(outfile_report, index = False, sep = '\t')

	# 9. Index bam files
    for spot_cell_type in ALL_Annotations:
    	outfile_bam="{}/{}.{}.bam".format(outdir, donor, spot_cell_type)
    	pysam.index(bam)
    
    print(f"Processing completed in {FILTER_dict['Total_time']} seconds.")


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------



# ============================================================================
# Function: initialize_parser
# ============================================================================


# Arguments

def initialize_parser():
	parser = argparse.ArgumentParser(description='Split alignment file into cell type specific BAMs')
	parser.add_argument('--bam', type=str, default=None, help='BAM file to be analysed (Sorted by coordinate)', required = True)
	parser.add_argument('--spot', type=str, default=None, help='File mapping spot barcodes to clusters information', required = False)
	parser.add_argument('--cell', type=str, default=None, help='File mapping cell barcodes to cell type information', required = False)
	parser.add_argument('--spot_cell', type=str, default=None, help='Metadata file mapping cell barcodes to spot barcodes', required = False)
	parser.add_argument('--id', type=str, default = 'Sample', help='Sample ID', required = False)
	parser.add_argument('--max_nM', type=int, default = None, help='Maximum number of mismatches permitted to consider reads for analysis. By default, this filter is switched off, although we recommed using --max_nM 5. If applied, this filter requires having the nM tag in the bam file. [Default: Switched off]', required = False)
	parser.add_argument('--max_NH', type=int, default = None, help='Maximum number of alignment hits permitted to consider reads for analysis. By default, this filter is switched off, although we recommend using --max_NH 1. This filter requires having the NH tag in the bam file. [Default: Switched off]', required = False)
	parser.add_argument('--min_MQ', type=int, default = 255, help='Minimum mapping quality required to consider reads for analysis. Set this value to 0 to switch this filter off. --min_MQ 255 is recommended for RNA data, and --min_MQ 30 for DNA data. [Default: 255]', required = False)
	parser.add_argument('--n_trim', type=int, default = 0, help='Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read [Default: 0]', required = False)
	parser.add_argument('--outdir', default = '.', help='Out directory', required = False)
	return (parser)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# ============================================================================
# Function: main
# ============================================================================


# Arguments

def main():
	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()
	bam = args.bam
	outdir = args.outdir
	spot_txt = args.spot
	cell_type_txt = args.cell
	spot_cell_map = args.spot_cell
	donor = args.id
	tissue = None
	max_NM = args.max_nM
	min_MAPQ = args.min_MQ
	n_trim = args.n_trim
	max_NH = args.max_NH

	# 2. Split bam file
	split_bam(bam, outdir, donor, tissue, max_NM, max_NH, min_MAPQ, n_trim, spot_txt, cell_type_txt, spot_cell_map)
    


#-------------
# Execute code
#-------------

if __name__ == '__main__':
	main()


