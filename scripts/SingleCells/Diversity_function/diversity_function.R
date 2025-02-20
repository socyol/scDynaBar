
# ==============================================================
#          🧬      DIVERSITY CALCULATION      🧬
# ==============================================================
# This script explains how the Diversity metric is calculated
# from barcode sequences in CRISPR barcoding experiments.
# 
# 🎯 Diversity is computed based on a pairwise alignment score 
# between the barcode sequence and the original reference spacer (gRNA).
# The higher the alignment score, the lower the diversity.
# To normalize this, we transform it by multiplying by -1.
#
# 🔬 Finally, we compute the **mean diversity per cell**.


# ==============================================================
# 📂 1. INPUT DATA
# ==============================================================
# - **barcode_sequences_10_Cells.csv**: sample of 10 cells with allele sequences
#   Columns: Cell_ID, Barcode_sequence, mean_Illumina_score
data <- read_csv("scDynaBar/scripts/SingleCells/Diversity_function/barcode_sequences_10_Cells.csv")

# 
# 📂 2. OUTPUT DATA
# ==============================================================
# - **cell_diversity_results.csv**: Table with per-cell diversity metrics
#   Columns: Cell_ID, Coverage, Diversity, Mean_length, Uncuts
output_path <- "scDynaBar/scripts/SingleCells/Diversity_function/"


#  Load required libraries
# ============================
library(Biostrings)
library(dplyr)
library(readr)

# ============================
# 📌 2. Define substitution matrix for alignment and the function
# ============================
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

alignment_function <- function(row) {
  sequence <- row[["Barcode_sequence"]]
  spacer <- DNAString("GGTATGCGGATGCAATCTCCGGGG")  # Reference gRNA spacer
  seq_to_align <- DNAStringSet(as.character(sequence))  
  
  if (!exists("seq_to_align", envir = environment())) {
    return(NULL) 
  } else {
    # Perform pairwise alignment (global-local)
    alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global-local", 
                                   substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
    
    if (!exists("alignment", envir = environment())) {
      return(NULL)
    } else {
      row$Alignment_score <- alignment@score  # Store alignment score
      row[["Length"]] <- nchar(sequence)  # Store barcode length
      row$Uncuts <- ifelse(as.character(alignment@subject) == as.character(spacer), "YES", "NO") # Identify uncut sequences
    }
  }
  return(row)
}

# ============================
# 📌 3. Read input data and apply alignment_function
# ============================
# Load example dataset: 100 barcode sequences (modify path if needed)
data <- read_csv("barcode_sequences_100.csv")
data<-r_2
result_list <- lapply(1:nrow(data), function(i) alignment_function(data[i, ]))
allele_table_features <- do.call(rbind, result_list)

# ============================
# 📌 4. Diversity function
# ============================
diversity_function <- function(data) {
  
  #  Remove rows with length constraints and NA values
  data <- subset(data, Length > 4 & Length < 35)
  data <- na.omit(data)
  
  #  Compute Diversity as the inverse of the alignment score
  data$Diversity <- -data$Alignment_score  # Inverting the alignment score
  data$Diversity <- (data$Diversity - min(data$Diversity)) / 
    (max(data$Diversity) - min(data$Diversity))  # Normalize
  
  #  Filter based on Illumina quality score
  data_filtered <- data[data$mean_Illumina_score >= 28, ]
  return(data_filtered)
}

data_filtered <- diversity_function(allele_table_features)

# ============================
# 📌 5. Compute mean Diversity per Cell
# ============================
process_cells_output <- function(data_filtered) {
  cells_output <- data_filtered %>%
    group_by(Cell_ID) %>%
    summarise(
      Coverage = n(),  # Number of reads per cell
      Diversity = mean(Diversity, na.rm = TRUE),  # Mean Diversity
      Mean_length = mean(Length, na.rm = TRUE),  # Mean barcode length
      Uncuts = sum(Uncuts == "YES") / n() * 100  # Percentage of uncut sequences
    )
  
  return(cells_output)
}

# Compute final table per Cell
cells_output <- process_cells_output(data_filtered)



# Save results
# ============================
write_csv(cells_output, file.path(output_path, "cell_diversity_results.csv"))


# 🎯 FINAL OUTPUT: "cell_diversity_results.csv" contains:
# - **Coverage**: Number of barcode sequences per cell
# - **Diversity**: Mean diversity (computed from alignment scores)
# - **Mean_length**: Average length of barcode sequences
# - **Uncuts**: Percentage of sequences that remain uncut


