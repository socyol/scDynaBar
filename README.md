#   🧬 scDynaBar 🧬

scDynaBar is an innovative approach that combines CRISPR-Cas9 dynamic barcoding with single-cell sequencing to record temporal cellular events. Over a 4-week period, genetic barcodes accumulate mutations, which are then sequenced together with the transcriptome of each single cell. This enables the creation of a time-ordered record of cellular events, providing a unique perspective on biological dynamics.

In this study, we applied scDynaBar to track the transition from a pluripotent state to a two-cell (2C)-like state in mouse embryonic stem cells (mESCs). Our results demonstrate the transient nature of the 2C-like state. Additionally, we show consistent mutation rates across different cell types in a mouse gastruloid model, underscoring the robustness and versatility of the system across diverse biological contexts.

<div style="display: flex; align-items: center;">
    <img src="https://github.com/user-attachments/assets/4baa7786-9729-4ab0-a8cd-b3fa9e6d4db1" alt="Image" width="150" style="margin-right: 10px;"/>
    <p style="font-size: 12px;">Scheme representation of the novel scDynaBar method. The
self-targeting guide RNA recognizes itself as it is in the spacer within the sequence of the
cassette inserted into each cell. Cas9 induces double-strand breaks that are repaired by
NHEJ. The Cassette or CRISPR-Cas9 barcode accumulates indels (red marks) over time.</p>
</div>



## Overview

This repository contains data and scripts for analyzing single-cell RNA-seq and Bulk RNA-seq experiments. The aim is to process barcode sequences and associated metadata, create Seurat objects, and generate visualizations for results presented in the associated paper.

## 📁 Repository Structure

### `data/`

- **`barcode_sequences/`**: This directory contains CSV files with barcode sequence data. Each file includes information on:
  - `cellID`: Identifier for each cell.
  - `barcode_sequence`: The specific barcode sequence (or allele).
  - `mean_illumina_score`: Quality score from Illumina sequencing for each single-cell experiment.

- **`metadatas/`**: This directory contains metadata CSV files that provide details on:
  - The number of reads/alleles (coverage) for each cell.
  - Allele features such as mean diversity, mean length, and percentage of original sequences (% original sequences).
  - For bulk experiments, this data is provided for each sample, including the system used (cas9 or BE3) and the spacer utilized (from a selection of 7).

- **`seurat_objects/`**: Contains the Seurat objects created from processed single-cell data, specifically those that have passed quality control.

### `scripts/`

#### Single-cells experiment analysis:
- **`sc_1_SeuratObject_analysis.R`**: This script processes the filtered_feature_bc_matrix to convert it into a Seurat object (matrix in GEOomnibus accession)

- **`sc_2_Barcode_sequences_analysis.R`**: This script processes the barcode sequence data and merges it with the Seurat objects to create the metadata (located in the `metadata` folder). This metadata is crucial for conducting analyses, viewing results, and generating figures.


#### Bulk analysis
- **`bulk_1_analysis.R`**: First step of the Bulk analysis. This script prepares the FASTQ data by organizing it into different folders to optimize and facilitate the subsequent analysis. This step ensures proper structuring and management of the data for the following phases.

- **`bulk_2_analysis.sh`**: Second step of the Bulk analysis. This Bash script automates the analysis process by creating jobs for each data file, which then execute the code provided in bulk_3_analysis.R. This parallelization helps to speed up the processing. At the end of this step, an output file similar to the barcode sequences for bulk data is generated, providing a comprehensive summary of the data processing.

- **`bulk_3_analysis.R`**:Third and final step of the Bulk analysis. This R script conducts a thorough analysis of the previously processed and organized data. It includes computations and data filtering.


#### PLOTS

- **`plots_1-bulk.R`**: Contains scripts for creating visualizations related to bulk data.

- **`plots_2-timecourse.R`**: Scripts for visualizations specific to time course analyses.

- **`plots_3-zscan4.R`**: Scripts for visualizations related to zscan4 experiment.
  
-  **`plots_4-gastruloids.R`**: Scripts for visualizations related to gastruloid data.

- **`settings.R`**: This script includes all necessary libraries and custom functions created for this project.


## 🛠️ How to run the Analysis
### Requirements ⚠️

- **R Version**: R version 4.3.2
  
1. Clone the repository 

```r
git clone https://github.com/yourusername/scDynaBar.git
cd scDynaBar
```
2. Install the required R packages:

```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
install.packages(c("Seurat", "dplyr", "patchwork", "Matrix", "ggplot2", 
                   "umap", "Rtsne", "gridExtra", "RColorBrewer", "stringr", 
                   "readxl", "openxlsx", "data.table", "SeuratObject", "psychTools"))

# Bioconductor packages
BiocManager::install(c("Biostrings", "ComplexHeatmap", "circlize"))

```
3. Load the necessary functions: Each script depends on helper functions and settings loaded from settings.R. Ensure that settings.R is sourced at the beginning of your analysis:
```r
source("scripts/settings.R")
```
4. Run the analysis scripts
5. Generate visualizations

## 📊 Raw data availability
The raw data will be made available through a GEO accession number on the Gene Expression Omnibus (GEO). The code will be updated once the raw data is accessible there, as it is too large to host directly in this repository.

## 📫 Contact
For questions or further assistance, please contact Yolanda Andrés-López at yalbmc@ibmb.csic.es.
