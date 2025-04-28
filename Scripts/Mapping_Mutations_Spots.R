library(Rsamtools)
library(dplyr)
library(tidyr)
library(GenomicRanges)

# ==============================================================================
#                             1. pVACseq data
# ==============================================================================

pvac_data <- read.csv("C:/Users/George_K/Desktop/Spatial_Statistics/test_output/MHC_Class_I/bam_GT.all_epitopes.aggregated.tsv", stringsAsFactors = FALSE, sep = '\t')

mutations <- pvac_data %>%
  
  select(ID, Gene, `AA.Change`, `Best.Peptide`, `IC50.MT`, Tier) %>%
  
  distinct() %>%
  
  filter(Tier == 'Pass', `IC50.MT` <= 50)
  
  

tmp = lapply(mutations$ID, function(x){
  
  x = strsplit(x, split = '-') |> unlist()
  names(x) = c('chr', 'Start', 'End', 'REF', 'ALT')
  
  return(x)
  
})

tmp = as.data.frame(do.call(rbind, tmp))

mutations = cbind.data.frame(tmp, mutations)

mut_granges <- makeGRangesFromDataFrame(mutations, keep.extra.columns = TRUE)


# ==============================================================================
#                             2. BAM file
# ==============================================================================

# Define parameters for scanning BAM
param <- ScanBamParam(
  what = c("pos", "qname", "seq", "qual", "cigar", "strand"),
  tag = c("CB"), # Cell barcode tag
  which = mut_granges # Only scan regions with mutations
)


mutation_mapping_function = function(bamfile, param){
  
  bam <- scanBam(file = bamfile, param = param)
  
  # Match mutations
  results <- data.frame()
  
  for (i in seq_along(bam)) {
    
    if (length(bam[[i]]$pos) > 0) {
      # Get mutation info
      mut_info <- mutations[i, ]
      
      # Get reads covering this mutation
      reads <- data.frame(
        spot_barcode = bam[[i]]$tag$CB,
        position = bam[[i]]$pos,
        cigar = bam[[i]]$cigar,
        sequence = bam[[i]]$seq,
        stringsAsFactors = FALSE
      )
      
      # Filter for reads that actually contain the mutation
      reads <- reads %>%
        mutate(
          read_pos = position,
          mut_pos_in_read = as.numeric(mut_info$Start) - read_pos + 1,
          contains_mutation = ifelse(
            mut_pos_in_read > 0 & mut_pos_in_read <= nchar(sequence),
            substr(sequence, mut_pos_in_read, mut_pos_in_read) == mut_info$ALT,
            FALSE
          )
        ) %>%
        filter(contains_mutation) %>%
        select(spot_barcode) %>%
        distinct() %>%
        mutate(
          chr = mut_info$chr,
          position = mut_info$Start,
          ref = mut_info$REF,
          alt = mut_info$ALT,
          gene = mut_info$Gene
        )
      
      results <- bind_rows(results, reads)
    }
  }
  
  return(results)
  
  
}

results = mutation_mapping_function(bamfile = 'C:/Users/George_K/Desktop/Spatial_Statistics/Data/Prostate_Cancer/Tumor_output.bam',
                                    param = param)

df_spots = read.csv("C:/Users/George_K/Desktop/Spatial_Statistics/Data/Prostate_Cancer/Spots_Barcodes.csv")

df_spots$Neoantigen_Status = ifelse(df_spots$Barcode %in% results$spot_barcode, 'Positive', 'Negative')


write.table(df_spots, 'Neoantigen_Status_pVACseq.csv', sep = ';', row.names = F)


# Count mutations per spot
spot_summary <- results %>%
  group_by(spot_barcode) %>%
  summarise(
    num_mutations = n(),
    mutated_genes = paste(unique(gene), collapse = ", "),
    mutation_details = paste(paste0(gene, ":", chr, ":", position, ":", ref, ">", alt), collapse = "; ")
  )

# Count spots per mutation
mutation_summary <- results %>%
  group_by(chr, position, ref, alt, gene) %>%
  summarise(
    num_spots = n(),
    spot_barcodes = paste(unique(spot_barcode), collapse = ", ")
  )




# # Save results
# write.csv(spot_summary, "visium_spots_with_mutations.csv", row.names = FALSE)
# write.csv(mutation_summary, "mutations_in_visium_spots.csv", row.names = FALSE)
# 
# # Print summary
# cat("Found", nrow(spot_summary), "spots containing", nrow(mutation_summary), "unique mutations\n")
# cat("Results saved to visium_spots_with_mutations.csv and mutations_in_visium_spots.csv\n")



