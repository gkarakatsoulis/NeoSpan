# NeoSpan
Neoantigen Prediction using Spatial Multi-Omic data - Pipeline

## ABOUT
This pipeline allows end-to-end analysis for neoantigen prediction using Spatial Multi-Omic data.
It incorporates spatial information along with pathologist region annotation (e.g. tumor vs normal) to identify, highlight and rank/evaluate potential neoantigens.

The pipeline can be divided into three primary categories:
* Mutation calling (SNV detection) using pathologist region annotation and spatial information
* Neoantigen prediction
* Ranking and evaluation of the detected potential neoantigens.

## Features
- Split/Filter the BAM file based on the regions annotated by the pathologist. It is strongly recommended to include normal (non-tumor) regions.
- Within each region, mutation calling (SNV detection) using spatial and/or cell-type information (Modified **SComatic** tool).
- Construct region specific feature-spot matrices and exclude those that appear in the normal region.
- Neoantigen prediction using the **pVACseq**.
- DGE analysis between mutated and unmutated spots, within the tumor region.
- Spatial statistics.

## Data requirements
- **BAM** file with a barcode tag (spot and/or cell).
- **Pathologist annotation (csv)**: Map each barcode to a specific region.
- **Spot clustering (csv)**: (Optional) To perform mutation calling within clusters of spots
- **Cell type annotation (csv)** (Optional) Only possible when we have single cell information apart from the spatial omics
- **Spot to cell mapping (csv)** (Optional) Only needed if there is both spot and cell-type information

## Installation
Create a Conda environment:
```bash
conda create -n vNS python=3.10
conda activate vNS
```

# Mutation calling (SNV detection) using pathologist region annotation, spatial and/or single cell information.

We show below how to run NeoSpan for the detection of SNVs using pathologist region annotation. It is a "two-stage" approach, which will be described bellow:

## Step 1: Split/Filter the BAM file based on the pathologist region annotation.
This step splits the BAM file into region-specific BAM files. To do so, the user provides a file mapping each spot barcode to a specific region. The script reads the unique regions that appear in the file, creates a new BAM file for each one separately, and writes it in an output directory provided by the user.

The step requires two data types as input:
* Aligned sequencing reads in BAM format for all spots analysed. The input BAM file must contain the spot barcode information in the tag “CB” (as reported by 10x Genomics).
* A csv file mapping each spot (named Barcode) to a specific region (named Region). Needs to be ";" delimitered. **Importanly, do not include doublets in this file**

**Execution**
```bash
Scripts/1_Split_Bam_Regions.py

```

## Step 2: Mutation calling within each region using the SComatic tool.
This step applies the *SComatic tool* (with slight modifications) to the region-specific BAM files created in Step 1. It considers spatial and/or cell types information.
To briefly describe it, it includes:
a) Splitting alignment file into spot clusters and/or cell-type-specific bam files
b) Collecting base count information
c) Merging base count matrices
d) Detection of somatic mutations

For a more thorough description, refer to the [SComatic README](https://github.com/cortes-ciriano-lab/SComatic/blob/main/README.md#detection-of-somatic-mutations-in-single-cell-data-sets-using-scomatic).

It requires the following data types as input:
* The region-specific BAM files created in step 1. The input files must contain the spot barcode information in the tag “CB”.
* (Optional) A csv file mapping each spot barcode to a specific cluster. The clustering analysis should have been conducted based on expression data.
* (Optional) A csv file mapping each cell barcode to a cell type. The annotation should have been conducted within each spot. Only available when both spatial and single-cell information is available.
* (Optional) A csv file mapping each cell barcode to a specific spot. Only needed if both spatial and single-cell information is available.

**Execution**
- Split the region-specific BAM files based on spot clusters and/or cell types.
```bash
Scripts/2_SplitBamSpotsCellTypes.py
```

- Construct a base count matrix for each spot cluster and/or cell type.
```bash
Scripts/3_BaseCellCounter.py
```

- Merge the count matrices derived above
```bash
Scripts/4_MergeBaseCellCounts.py
```


- Apply filters and Beta binomial tests to discount sites affected by recurrent technical artefacts as somatic mutations.

```bash
Scripts/5a_BetaBinEstimation.py # Estimate the Beta Binomial parameters

Scripts/5b_BaseCellCalling.step1.py
```

- Apply additional filters based on external datasets (RNA editing and Panel of Normals).
```bash
Scripts/5c_BaseCellCalling.step2.py

```

# Neoantigen prediction

We show below how to run NeoSpan for the Neoantigen Prediction. The procedure includes the following steps:
1) Define the normal (non-tumor) region. For this reason, it is strongly recommended to include normal regions.
2) For the tumor region(s), remove the SNVs that were also detected in the normal region. This step is critical since by definition the neoantigens are antigens generated by tumor cells. Ignore it if no normal regions exist in the dataset.
2) Construct the feature-spot matrices for each region separately.
3) Apply the pVACseq 
