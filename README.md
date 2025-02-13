# NeoSpan
Neoantigen Prediction using Spatial Multi-Omic data - Pipeline

## ABOUT
This pipeline allows end-to-end analysis for neoantigen prediction using Spatial Multi-Omic data.
It incorporates spatial information along with pathologist region annotation (tumor vs normal) to identify, highlight and rank/evaluate potential neoantigens.

The pipeline can be divided into three primary categories:
* Mutation calling (SNV detection) using pathologist region annotation and spatial information
* Neoantigen prediction
* Ranking and evaluation of the detected potential neoantigens.

## Features
- Split/Filter the BAM file based on the Tumor and Normal regions.
- Within each region, mutation calling (SNV detection) using spatial and/or cell-type information (Modified **SComatic** tool).
- Construct region specific feature-spot matrices and exclude those that appear in the normal region.
- Neoantigen prediction using the pVACseq.
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

