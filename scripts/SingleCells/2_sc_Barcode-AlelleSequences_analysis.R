##############################################
###   Barcode-sequences analysis  🧬 📊    ###
##############################################

#   This script is to analyze the scDynaBar Alleles sequences data frames from the 
#   single cells experiments: timecourse, gastruloids and zscan4.
#
#   We process the table of alleles for each cell and apply different
#   filters to conduct the analysis and study. 

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R") # script with ALL the functions

#_______________________________________
# 2. Load Barcode_sequence file 
# ======================================
input_bs <- "data/barcode_sequences/" # folder of the allele sequences
barcode_sequence <- read.csv(file.path(input_bs, "BARCODE_SEQ_FILE")) # ⚠️ select the file depending on the experime️nt  ⚠️

output_path <- "results/" # save the result after applying the function and filters
if (!dir.exists(output_path)) dir.create(output_path)


#_______________________________________
# 3. Analyze and filter the sequences table
# ======================================

# Apply the function to perform Pairwise alignment with the original gRNA 
result_list <- lapply(1:nrow(barcode_sequence), function(i) read_table_function_modified(barcode_sequence[i, ])) 
read_table <- do.call(rbind, result_list)

saveRDS(read_table, file = file.path(output_path, "NAME_READTABLE")) # you decide the filename
#   >>> this dataframe is useful to represent the supplementary figure 5 about
#       illumina score threshold for timecourse and zscan4 experiment 

# Filter the barcode sequences or alleles 
filtered_data <- process_data(read_table)

# Create a table by cell_ID with coverage, mean diversity, mean length and %Uncuts
cell_ID_table <- process_cells_output(filtered_data)
cell_ID_table <- cell_ID_table[cell_ID_table$Coverage > 3,]


#_______________________________________
# 4. Merge with seurat object
# ======================================

#         >>>>>>    We want to merge the info from the Barcode/Alleles table
#                   with the Seurat Objects to create metadatas and perform
#                   the plots for the paper 

#                   

# Load the seurat object already filtered (depending on the experiment)
input_so <- "data/seurat_objects/"
sce <- readRDS(file.path(input_so, "SEURAT_OBJECT.rds")) # select the seurat object

# select the cells in common with the cell_ID_table and filter
colnames_to_keep <- intersect(colnames(sce), cell_ID_table$Cell_ID)
sce_filtered <- sce[, colnames_to_keep]
sce_filtered@meta.data$Cell_ID <- rownames(sce_filtered@meta.data)

# Merge metadata by Cell_ID
sce_filtered@meta.data <- merge(sce_filtered@meta.data, cell_ID_table, by = "Cell_ID")
rownames(sce_filtered@meta.data) <- sce_filtered@meta.data$Barcode 

# From here we can extract the metadata from the merged_seuratObject
# or built the metadata (see metadata/ folder)
# to perform the plots and figures in the paper 
#  input_metadata <- "data/metadatas/"