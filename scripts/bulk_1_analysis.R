###############################
#       BULK ANALYSIS 1
###############################
library(data.table)
library(ShortRead) # BiocManager::install("ShortRead")


## 1) Load the table and the directories 
# ---------------------------------------------------
fastq_dir <- "~YOUR_PATH/FASTQ_FOLDER" 
output_dir <- "~YOUR_PATH/Processed_CSVs/"  # directory to save csv files
# create the directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}



## 2) Define a function to process the FASTQ File
# ---------------------------------------------------
process_fastq_to_df <- function(fastq_file) {
  
  fq <- readFastq(fastq_file)
  sequences <- as.character(sread(fq)) # extract nt sequence
  quality_scores <- as.character(as(quality(fq), "PhredQuality"))
  
  df <- data.frame(
    Read = seq_along(sequences),        
    Sequence = sequences,            
    Quality_sequence = quality_scores,     
    stringsAsFactors = FALSE          
  )
  
  return(df)
}

## example:
# fasqc <- "ane8227_AACAACC_A06_L001_R1.fastq.gz"
# fastq_df <- process_fastq_to_df(fasqc)
# View(fastq_df)

## 3) List of all the FASTQ files and APPLY the Function
# ---------------------------------------------------
fastq_files <- list.files(path = fastq_dir, full.names = TRUE)

count = 0
for (fastq_file in fastq_files) {
  processed_df <- process_fastq_to_df(fastq_file)
  processed_df$File <- basename(fastq_file)
  exp_name <- paste0(
    gsub(".*_([A-Za-z][0-9]{1,2})_.*", "\\1", basename(fastq_file)), 
    "_",
    gsub(".*lane([0-9]+).*", "\\1", basename(fastq_file))
  )
  
  file_name <- paste0(exp_name, ".csv")
  file_path <- file.path(output_dir, file_name)
  
  write.csv(processed_df, file = file_path, row.names = FALSE)
  count = count+1
  
  cat(count, "Archivo CSV guardado en:", file_path, "\n")
}


## 4) Create folders for each gRNA
# ---------------------------------------------------
csv_dir <- output_dir
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)

valid_csv_files <- bulk$V1
length(valid_csv_files)

output_base_dir <- "~/Documents/TFM/Bulk_analysis/gRNA_Folders/" # directory to save all the files
if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir)
}

for (i in 1:nrow(bulk)) {
  file_name <- bulk$V1[i] # name of the file
  matching_file <- grep(file_name, csv_files, value = TRUE)
  
  if (length(matching_file) > 0) {
    gRNA <- bulk$gRNA[i] # obtain the gRNA of the file (if the file is found)
    # create the folder for that gRNA
    output_dir <- file.path(output_base_dir, gRNA)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    file.copy(matching_file, file.path(output_dir, basename(matching_file))) # move the file to the correspondent folder (grna)
    
    cat("File", matching_file, "copied in folder", output_dir, "\n")
  } else {
    cat("File", file_name, "not found in csv_files\n")
  }
}
