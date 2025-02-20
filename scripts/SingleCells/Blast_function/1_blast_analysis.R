
# ==============================================================
#           OFF-TARGET ANALYSIS OF gRNAs WITH BLAST
# ==============================================================
# This script processes the gRNA sequences from the Day 4 and Day 10 
# CRISPR barcoding experiment to generate a FASTA file for BLAST 
# and a CSV file summarizing the gRNA frequencies.
#
# 📌 Steps:
# 1️⃣ Load the data from a CSV file containing barcode allele information.
# 2️⃣ Compute the frequency of each unique gRNA sequence.
# 3️⃣ Save the sequences in a FASTA format for BLAST analysis.
# 4️⃣ Generate a CSV file summarizing the gRNA sequences and their frequencies.


# LOAD REQUIRED LIBRARIES
# ==============================================================
library(dplyr)
library(readr)

# ==============================================================
#   FUNCTION TO GENERATE gRNA FASTA & CSV FILES
# ==============================================================
generate_gRNA_fasta_and_csv <- function(input_csv, output_fasta, output_csv) {
  
  reads <- read_csv(input_csv)
  
  match_counts <- table(reads$Match)
  total_reads <- sum(match_counts)  
  match_counts <- sort(match_counts, decreasing = TRUE)
  sequences <- names(match_counts)
  
  file_conn <- file(output_fasta, open = "w")
  
  csv_data <- data.frame(gRNA = character(), Spacer = character(), Frequency = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(sequences)) {
    freq <- match_counts[sequences[i]] / total_reads 
    
    writeLines(paste0(">gRNA", i, ",", round(freq * 100, 4), "%"), file_conn)
    writeLines(sequences[i], file_conn)
    
    csv_data <- rbind(csv_data, data.frame(
      gRNA = paste0("gRNA", i),
      Spacer = sequences[i],
      Frequency = round(freq * 100, 4) 
    ))
  }
  
  close(file_conn)
  
  write_csv(csv_data, output_csv)
  
  # Print success messages
  cat("✅ FASTA file created:", output_fasta, "\n")
  cat("✅ CSV file created:", output_csv, "\n")
}


# ==============================================================
#  INPUT DATA & FILE PATHS and RUN THE FUNCTION TO GENERATE FASTA & CSV
# ==============================================================
# The input file should contain barcode-allele match information.
# Ensure the file has a "Match" column representing gRNA sequences.

input_csv <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/barcodes-alleles_cuts_day4-10.csv"
output_fasta <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/gRNA_sequences_with_frequencies_DAY4-10.fasta"
output_csv <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/gRNA_sequences_with_frequencies_DAY4-10.csv"

generate_gRNA_fasta_and_csv(input_csv, output_fasta, output_csv)

reads <- read_csv(output_csv)
cat("📊 Number of gRNA sequences processed:", nrow(reads), "\n")

# ==============================================================
# 🔗 NEXT STEP: 2_BLAST_ANALYSIS.sh
# ==============================================================
# The generated FASTA file can now be used for BLAST analysis.
# The second step of BLAST processing will be handled separately.





















################################
##  1.  FUNCTION WITH ALL gRNA
# ================================================
generate_gRNA_fasta_and_csv <- function(input_csv, output_fasta, output_csv) {
  reads <- read.csv(input_csv)
  match_counts <- table(reads$Match)
  total_reads <- sum(match_counts)  
  
  
  match_counts <- sort(match_counts, decreasing = TRUE)
  sequences <- names(match_counts)
  
  
  file_conn <- file(output_fasta, open = "w")
  
  csv_data <- data.frame(gRNA = character(), Spacer = character(), Frequency = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(sequences)) {
    freq <- match_counts[sequences[i]] / total_reads 
    
    writeLines(paste0(">gRNA", i, ",", round(freq * 100, 4), "%"), file_conn)
    writeLines(sequences[i], file_conn)
    
    csv_data <- rbind(csv_data, data.frame(
      gRNA = paste0("gRNA", i),
      Spacer = sequences[i],
      Frequency = round(freq * 100, 4) 
    ))
  }
  
  close(file_conn)
  
  write.csv(csv_data, output_csv, row.names = FALSE)
  
  cat("FASTA file created:", output_fasta, "\n")
  cat("CSV file created:", output_csv, "\n")
}


##  2.  Run the function
# ............................................
input_csv <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/barcodes-alleles_cuts_day4-10.csv"
output_fasta <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/gRNA_sequences_with_frequencies_DAY4-10.fasta"
output_csv <- "Documents/TFM/scDynaBar/scripts/SingleCells/Blast_function/gRNA_sequences_with_frequencies_DAY4-10.csv"

# Generar el archivo FASTA y el archivo CSV
generate_gRNA_fasta_and_csv(input_csv, output_fasta, output_csv)


# reads <- read.csv("Documents/TFM/paper_REVIEWER-COMENTS/gRNA_sequences_with_frequencies_DAY4-10_30-01-25.csv")
dim(reads) # 1480

## 4. meterlo en el cluster y correr la función de Blast
# ..........................................................
#>>>> comento los comandos que he usado en bash;;
#
# scp -r gRNA_sequences_with_frequencies_DAY4-10_30-01-25.fasta yalbmc@10.5.60.113:/embcri-data/ihhbmc/yalbmc/barcodes/blast_paper

# he usado job_blast_30-01-25.sh