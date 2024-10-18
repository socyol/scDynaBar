# scDynaBar

Project Description

scDynaBar is an innovative approach that combines CRISPR-Cas9 dynamic barcoding with single-cell sequencing to record temporal cellular events. Over a 4-week period, genetic barcodes accumulate mutations, which are then sequenced together with the transcriptome of each single cell. This enables the creation of a time-ordered record of cellular events, providing a unique perspective on biological dynamics.

In this study, we applied scDynaBar to track the transition from a pluripotent state to a two-cell (2C)-like state in mouse embryonic stem cells (mESCs). Our results demonstrate the transient nature of the 2C-like state. Additionally, we show consistent mutation rates across different cell types in a mouse gastruloid model, underscoring the robustness and versatility of the system across diverse biological contexts.

## Overview

This repository contains data and scripts for analyzing single-cell experiments. The aim is to process barcode sequences and associated metadata, create Seurat objects, and generate visualizations for results presented in the associated paper.

## Repository Structure

### `data/`

- **`barcode_sequences/`**: This directory contains CSV files with barcode sequence data. Each file includes information on:
  - `cellID`: Identifier for each cell.
  - `barcode_sequence`: The specific barcode sequence (or allele).
  - `quality_score`: Quality score from Illumina sequencing for each single-cell experiment.

- **`metadatas/`**: This directory contains metadata CSV files that provide details on:
  - The number of reads/alleles (coverage) for each cell.
  - Allele features such as mean diversity, mean length, and percentage of original sequences (% original sequences).
  - For bulk experiments, this data is provided for each sample, including the system used (cas9 or BE3) and the spacer utilized (from a selection of 7).

- **`seurat_objects/`**: Contains the Seurat objects created from processed single-cell data, specifically those that have passed quality control.

### `scripts/`

- **`1_Barcode_sequences_analysis.R`**: This script processes the barcode sequence data and merges it with the Seurat objects to create the metadata (located in the `metadata` folder). This metadata is crucial for conducting analyses, viewing results, and generating figures.

- **`2_plots-bulk.R`**: Contains scripts for creating visualizations related to bulk data.

- **`3_plots-timecourse.R`**: Scripts for visualizations specific to time course analyses.

- **`4_plots-gastruloids.R`**: Scripts for visualizations related to gastruloid data.

- **`settings.R`**: This script includes all necessary libraries and custom functions created for this project.

## Requirements

- **R Version**: R version 4.3.2

### Required Packages

To run the scripts, ensure you have the following R packages installed:
