
#############################
#         FUNCTIONS
#############################

# Settings - utils
# Scripts with ALL the functions:

# Libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork)) 
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(psychTools))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))



#        ====================================
#           Barcode Sequences ANALYSIS  STEP   
#        ====================================

# ....................................
#     READS/Allele - TABLE function
# ....................................
# function reads table with the diversity and alignment in TimeCourse (Day 0, Day 4, Day 10), Zscan4c and Gastruloids (day 6) experiments
#       (in case of Bulk data we use different spacers)



sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

read_table_function_modified <- function(row) {
  sequence <- row[["Barcode_sequence"]]
  spacer <- DNAString("GGTATGCGGATGCAATCTCCGGGG")   # spacer
  seq_to_align <- DNAStringSet(as.character(sequence))  
  if (!exists("seq_to_align", envir = environment())) {
    return(NULL) 
  } else {
    
    # Perform a pairwise alignment with the barcode allele and the original spacer (gRNA)
    alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global-local", substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
    
    if (!exists("alignment", envir = environment())) {
      return(NULL)
    } else {
      row$Alignment_score <- alignment@score # save the alignment score to compute diversity feature
      row[["Length"]] <- nchar(sequence) # compute length of the barcode sequence
      row$Uncuts <- ifelse(as.character(alignment@subject) == as.character(spacer), "YES", "NO") # to later compute % of uncuts
    }
  }
  
  return(row)
}
#  Feature calculation
# =========================
# result_list <- lapply(1:nrow(allele_table), function(i) read_table_function_modified(allele_table[i, ]))
# allele_table_features <- do.call(rbind, result_list)


# ....................................
#    Processing/Filtering FUNCTION
# ....................................

#       Function for applying filters and compute mean of the mutation rate (Diversity) in the alleles table:

process_data <- function(data, output_path = "") {
  
  # 1) Remove rows with Length < 4 or > 35 and remove NA values
  data <- subset(data, Length > 4 & Length < 35)
  data <- na.omit(data)
  
  # 2) Add a "Diversity" column based on the negative of "Alignment_score"
  data$Diversity <- -data$Alignment_score
  data$Diversity <- (data$Diversity - min(data$Diversity)) / 
    (max(data$Diversity) - min(data$Diversity))
  
  # 3) Filter rows with Illumina_score >= 28
  data_filtered <- data[data$mean_Illumina_score >= 28, ]
  
  # Return the filtered data
  return(data_filtered)
}

# filtered_data <- process_data(reads2_filtered)


# .....................................................
#    # Create TABLE of alleles by barcode FUNCTION
# .....................................................

#       Function to create a table by CellID with features such as coverage, mean diversity, %Uncuts and Mean length:

process_cells_output <- function(data_filtered, output_path = "YOUR_PATH") {
  cells_output <- data_filtered %>%
    group_by(Cell_ID) %>%
    summarise(
      Coverage = n(),
      Diversity = mean(Diversity, na.rm = TRUE),
      Mean_length = mean(Length, na.rm = TRUE),
      Uncuts = sum(Uncuts == "YES") / n() * 100
    )
  
  #write.csv(cells_output, file = output_path, row.names = FALSE)
  return(cells_output)
}

# cells_output <- process_cells_output(filtered_data, "YOUR_PATH")

# Apply a filter of coverage > 3
# ==================================
# cells_output <- cells_output[cells_output$Coverage > 3,]






#######################################################################################
#  functions for each fQC filters (mt and nFeatures) to do the plots for the paper:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FilterAndPlotMT <- function(sce,  
                            max_percent_mt = 7.5, 
                            title = "Gastruloids 3", 
                            colors = c("#FFB6C1","skyblue3")) {
  
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
  sce_filtered_mt <- subset(sce, subset = percent.mt < max_percent_mt)
  sce@meta.data$filtered <- ifelse(rownames(sce@meta.data) %in% rownames(sce_filtered_mt@meta.data), "Retained", "Filtered")
  plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "filtered", cols = colors)
  p1 <- plot1 + 
    geom_hline(yintercept = max(sce_filtered_mt@meta.data$percent.mt), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Retained" = colors[2], "Filtered" = colors[1])) +
    ggtitle(title)
  theme(legend.position = "none")  # Ocultar la leyenda
  
  return(list(filtered_data = sce_filtered_mt, plot1 = p1))
}

FilterAndPlotRNA <- function(sce,  
                             min_counts = 3000, 
                             title = "Gastruloids 3", 
                             colors = c("#FFB6C1","skyblue3")) {
  sce_filtered_rna <- subset(sce, subset = nFeature_RNA > min_counts)
  sce@meta.data$filtered <- ifelse(rownames(sce@meta.data) %in% rownames(sce_filtered_rna@meta.data), "Retained", "Filtered")
  plot2 <- FeatureScatter(sce, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = "filtered", cols = colors)
  p2 <- plot2 + 
    geom_vline(xintercept = min(sce_filtered_rna@meta.data$nFeature_RNA), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Retained" = colors[2], "Filtered" = colors[1])) +
    ggtitle(title)
  theme(legend.position = "none")  # Ocultar la leyenda
  
  return(list(filtered_data = sce_filtered_rna, plot2 = p2))
}




#######################################################################################

##      🚨     The functions bellow is for RAW DATA  (GEO omnibus)      🚨


#        ====================================
#                  Preprocessing STEP
#        ====================================

# The functions bellow are for the 🚨 raw data (seurat object unprocessed) 🚨

# -----------------------------
#   Create Seurat Object 
# -----------------------------
sce.data <- Read10X(data.dir = "YOURPATH/filtered_feature_bc_matrix/")
colnames(sce.data) <- gsub("-1", "", colnames(sce.data))
sce <- CreateSeuratObject(counts = sce.data, project = "Single-cells", min.cells = 3, min.features = 200)

# FROM THE SCE CREATED TO APPLY ALL THE FUNCTIONS FOR QC

# ....................
#     QC function
# ....................
# This function is to perform the quality control and get the cells after the filter
QualityControl <- function(sce, max_percent_mt = 7.5, min_counts = 3000) {
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
  sce_filtered <- subset(sce, subset = 
                           percent.mt < max_percent_mt & 
                           nFeature_RNA > min_counts)
  n_cells <- ncol(sce)  
  n_filt_cells <- ncol(sce_filtered) 
  perc_filtered <- (n_filt_cells / n_cells) * 100
  cat("\n", "Number of cells before filtering =", n_cells, "\n", 
      "Number of cells after filtering =", n_filt_cells, "\n",
      "Percentage of filtered cells:", perc_filtered, "%\n",
      "Number of genes:", dim(sce_filtered), "\n")
  return(sce_filtered)
}

# sce_filtered <- QualityControl(sce)


# Function to extract only the sce filtered
metadata_mit <- function(sce, max_percent_mt = 7.5, min_counts = 3000) {
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
  
  return(sce)
}

# Function QC with plots 
FilterAndPlot <- function(sce,  
                          max_percent_mt = 7.5, 
                          min_counts = 3000, 
                          title = "Gastruloids 3", 
                          colors = c("#FFB6C1","skyblue3")) {
  
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
  sce_filtered <- subset(sce, subset = 
                           percent.mt < max_percent_mt & 
                           nFeature_RNA > min_counts)
  
  sce_filtered <- subset(sce, subset = 
                           percent.mt < max_percent_mt & 
                           nFeature_RNA > min_counts)
  
  sce@meta.data$filtered <- ifelse(rownames(sce@meta.data) %in% rownames(sce_filtered@meta.data), "Retained", "Filtered")
  
  plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "filtered", cols = colors)
  
  p1 <- plot1 + 
    geom_hline(yintercept = max(sce_filtered@meta.data$percent.mt), linetype = "dashed", color = "red") +
    # geom_vline(xintercept = min(sce_filtered@meta.data$nFeature_RNA), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Retained" = colors[2], "Filtered" = colors[1])) +
    ggtitle(title)
  
  plot2 <- FeatureScatter(sce, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = "filtered", cols = c("skyblue3", "#FFB6C1"))
  p2 <- plot2 + 
    # geom_hline(yintercept = max(sce_filtered@meta.data$percent.mt), linetype = "dashed", color = "red") +
    geom_vline(xintercept = min(sce_filtered@meta.data$nFeature_RNA), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Retained" = colors[2], "Filtered" = colors[1])) +
    ggtitle(title)
  
  return(list(filtered_data = sce_filtered, plot1 = p1, plot2 = p2))
}


# ..................................
# Normalize, Scale, and PCA function
# ..................................
NormalizeScalePCA <- function(sce, nfeatures = 2000, scale_factor = 10000) {
  
  # Normalization and Variable Feature Selection
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = scale_factor)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = nfeatures)
  
  # Get genes for scaling (use all genes available)
  genes <- rownames(sce@assays$RNA)
  
  # Scaling and PCA
  sce <- ScaleData(sce, features = genes)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  
  # Visualize PCA elbow plot
  ElbowPlot(sce)
  
  return(sce)
}

# Example of usage
# sce_combined <- NormalizeScalePCA(sce_filtered)


# ...........................
#     GFP POSITIVES
# ...........................
filter_gfp_positive_cells <- function(sce_norm) {  
  genes <- rownames(sce_norm@assays$RNA)  
  genes_gfp <- genes[grepl("^GFP", genes)]
  
  # Get the assay data for GFP-related genes
  genes_values <- GetAssayData(sce_norm, layer = "counts")[genes_gfp, ]  # or put "data"
  sce_norm@meta.data$gfp <- genes_values  
  total_cells <- ncol(sce_norm)  
  sce_gfp_positive <- subset(sce_norm, subset = gfp > 0)  
  gfp_positive_cells <- ncol(sce_gfp_positive)  
  percentage_gfp_positive <- (gfp_positive_cells / total_cells) * 100  
  cat("Percentage of GFP-positive cells:", percentage_gfp_positive, "%\n")  
  return(sce_gfp_positive)
}

# result <- filter_gfp_positive_cells(sce_norm)


# ---------------------------------------------------------------------



#        ====================================
#         Barcode Sequences ANALYSIS  STEP    >>>>> 🚨 FOR RAW DATA (GEO OMNIBUS) 🚨
#        ====================================
# The functions bellow are for the raw data of the allele sequences (GEOomnibus repository :))

# ....................................................
#    Filter Cassette/reads/Allele table function
# ....................................................
# Function to filter the raw output (the .txt file with the cassettes and quality of ilumina)
filter_output_amplicon <- function(file_path, gfp_positive_barcodes, sfl1 = "ACACC", sfl2 = "TTAGAG") {
  
  # 1. Load the output file
  output <- read.table(file_path, header = FALSE, sep = "\t")
  colnames(output)[colnames(output) == "V1"] <- "Barcode"
  colnames(output)[colnames(output) == "V2"] <- "Cassette"
  colnames(output)[colnames(output) == "V3"] <- "Quality"
  
  # 2. Apply filters
  output_filtered <- output[grepl(paste0(sfl1, ".*", sfl2), output$Cassette), ]
  output_filtered <- output_filtered[!grepl("N", output_filtered$Cassette), ]
  g <- output_filtered[output_filtered$Barcode %in% gfp_positive_barcodes, ]
  
  cat(" Number of cells in the original:", dim(output)[1], "\n",
      "Number of cells in the filtered:", dim(g)[1])
  return(g)
}

#g <- filter_output_amplicon(
#  file_path = "Documents/gastruloids/3-gastruloids_output_amplicon_gastruloid_3-txt_2024-07-31_1509/Output_amplicon_Gastruloid_3.txt", 
#  gfp_positive_barcodes = colnames(sce_gfp_positive)
#)


# ...........................
#     READS-TABLE function
# ...........................
# function reads table with the diversity and alignment in TimeCourse (Day 0, Day 4, Day 10), Zscan4c and Gastruloids (day 6) experiments
#       (in case of Bulk data we use different spacers)

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
read_table_function_modified <- function(row) {
  
  # 1- Extract the QUALITY part of the cassette (mean quality, quality cas)
  posicion1<- row[["POSICION_SFL1"]]
  posicion2<- row[["POSICION_SFL2"]]
  quality_text <- row[["Quality"]]
  
  substring_quality <- substr(quality_text, posicion1, posicion2) # + attr(posicion2, "match.length") - 1)
  row[["CAS_QUALITY"]] <- substring_quality
  
  quality_score <- row[["CAS_QUALITY"]]
  quality_scores_ascii <- lapply(quality_score, charToRaw)
  mean_scores <- sapply(quality_scores_ascii, function(x) mean(as.integer(x)))
  row[["Quality_Scores"]] <- mean_scores
  
  # 2- Extract the barcode cassettes in another column to perform alignment
  match_column <- row[["Match"]]
  row[["Length"]] <- nchar(match_column)
  
  # 3- GLOBAL ALIGNMENT
  spacer <- DNAString("GGTATGCGGATGCAATCTCCGGGG")  
  seq_to_align <- DNAStringSet(as.character(match_column))  
  
  if (!exists("seq_to_align", envir = environment())) {
    return(NULL)  # Retorna NULL para eliminar la row
  } else {
    
    # sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
    alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global-local", substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
    
    if (!exists("alignment", envir = environment())) {
      return(NULL)
    } else {
      row$Alineamiento_Score <- alignment@score  
      row$Original <- ifelse(as.character(alignment@subject) == as.character(spacer), "YES", "NO")   
    }
  }
  
  return(row)
}


#  Feature calculation
# =========================
# result_list <- lapply(1:nrow(subset_umi), function(i) read_table_function_modified(subset_umi[i, ]))
# read_table <- do.call(rbind, result_list)


#  Processing
# =======================
process_data <- function(data, output_path = "") {
  
  # 1) Remove rows with Length < 4 or > 35 and remove NA values
  data <- subset(data, Length > 4 & Length < 35)
  data <- na.omit(data)
  
  # 2) Add a "Diversity" column based on the negative of "Alineamiento_Score"
  data$Diversity <- -data$Alineamiento_Score
  data$Diversity <- (data$Diversity - min(data$Diversity)) / 
    (max(data$Diversity) - min(data$Diversity))
  
  data$Quality_Scores <- data$Quality_Scores - 33
  
  # 3) Filter rows with Quality_Scores >= 55
  data_filtered <- data[data$Quality_Scores >= 28, ]
  
  # Save the filtered data to a CSV file
  write.csv(data_filtered, file = output_path, row.names = FALSE)
  
  # Return the filtered data
  return(data_filtered)
}

# filtered_data <- process_data(reads2_filtered, "YOUR_PATH")





