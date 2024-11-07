##############################################
###      Seurat Object analysis        ###
##############################################

#   This script is to analyze and create the seurat object of the single cells
#   experiments (timecourse, gastruloids and zscan4) from the matrix published in GEOmnibus 
#   

#   We do QC, normalization, scalling, GFP selection and then this result will be aligned 
#   with the barcode sequences (alleles) analysis to conduct our study

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

#_______________________________________
# 2. Load filtered_feature_bc_matrix
# ======================================
# -----------------------------
#  Seurat Object  >>> QC
# -----------------------------
sce.data <- Read10X(data.dir = "YOUR_PATH/filtered_feature_bc_matrix/")
colnames(sce.data) <- gsub("-1", "", colnames(sce.data))
sce <- CreateSeuratObject(counts = sce.data, project = "Single-cells", min.cells = 3, min.features = 200)


sce_filtered <- QualityControl(sce)
sce_norm <- NormalizeScalePCA(sce_filtered)
sce_gfp_pos <- filter_gfp_positive_cells(sce_norm)

# TSNE or UMAP will be done depending on parameters in each experiment

#_______________________________________