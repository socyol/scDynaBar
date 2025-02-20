###############################
#       BULK ANALYSIS 3
###############################

# Function to analyse each Allele:

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

read_table_function <- function(row, spacer) {
  sfl1.1 <- "AACAC"
  sfl2.1 <- "TTAGAG"
  
  cassette <- row[["Sequence"]]
  quality_text <- row[["Quality_sequence"]]
  
  if (is.null(cassette) || cassette == "" || is.na(cassette)) {
    return(NULL)  # if there is not cassette --> NULL
  }
  
  positions_sfl1 <- gregexpr(sfl1.1, cassette)[[1]]
  positions_sfl2 <- gregexpr(sfl2.1, cassette)[[1]]
  
  if (length(positions_sfl1) > 1) {
    positions_sfl1 <- positions_sfl1[1]
  }
  
  if (length(positions_sfl2) > 1) {
    positions_sfl2 <- positions_sfl2[1]
  }
  
  if (positions_sfl1 == -1 || positions_sfl2 == -1) {
    return(NULL) 
  }
  
  row[["POSICION_SFL1"]] <- positions_sfl1
  row[["POSICION_SFL2"]] <- positions_sfl2
  row[["SFL1"]] <- 1
  row[["SFL2"]] <- 1
  
  if (is.na(positions_sfl1) || is.na(positions_sfl2) || positions_sfl1 < 1 || positions_sfl2 < 1) {
    return(NULL)  
  }
  
  quality_length <- nchar(quality_text)
  
  match_length <- attr(positions_sfl2, "match.length")
  if (is.null(match_length)) {
    match_length <- 0 
  }
  
  if (positions_sfl1 > quality_length || (positions_sfl2 + match_length - 1) > quality_length) {
    return(NULL)  
  }
  
  
  substring_quality <- substr(quality_text, positions_sfl1, positions_sfl2 + match_length - 1)
  row[["CAS_QUALITY"]] <- substring_quality
  
  quality_scores_ascii <- lapply(row[["CAS_QUALITY"]], charToRaw)
  mean_scores <- sapply(quality_scores_ascii, function(x) mean(as.integer(x)))
  row[["QS_Illumina"]] <- mean_scores - 33
  
  row[["Allele"]] <- gsub(paste0(".*(", sfl1.1, ".*", sfl2.1, ").*"), "\\1", cassette)
  row[["Allele"]] <- gsub(paste0(".*(?:", sfl1.1, ")(.*?)(?:", sfl2.1, ").*"), "\\1", cassette)
  row[["Length"]] <- nchar(row[["Allele"]])
  
  # Global-local alignment 
  seq_to_align <- DNAStringSet(as.character(row[["Allele"]]))
  alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global-local", substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
  row$Alignment_score <- alignment@score
  row$Uncuts <- ifelse(as.character(alignment@subject) == as.character(spacer), "YES", "NO")
  subject <- toString(alignedSubject(alignment))
  pattern <- toString(alignedPattern(alignment))
  
  row$insertions <- sum(strsplit(pattern, "")[[1]] == "-")
  row$deletions <- sum(strsplit(subject, "")[[1]] == "-")
  row$substitutions <- nmismatch(alignment)
  
  return(row)
}



# Process the result
process_file <- function(input_file, output_dir, spacer) {
  
  d <- read.csv(input_file)
  sfl1 <- "AACAC"
  sfl2 <- "TTAGAG"
  d_filtered <- d[grepl(paste0(sfl1, ".*", sfl2), d$Sequence), ]
  d_filtered <- d_filtered[!grepl("N", d_filtered$Sequence), ]
  resultado_lista <- lapply(1:nrow(d_filtered), function(i) read_table_function(d_filtered[i, ], spacer))
  resultado <- do.call(rbind, resultado_lista)
  resultado$Diversity <- -resultado$Alignment_score
  
  if (!all(is.na(resultado$Diversity))) {
    resultado$Diversity <- (resultado$Diversity - min(resultado$Diversity, na.rm = TRUE)) / 
      (max(resultado$Diversity, na.rm = TRUE) - min(resultado$Diversity, na.rm = TRUE))
  }
  
  output_file <- paste0(output_dir, "/", basename(input_file), "_table_alleles.csv")
  write.csv(resultado, output_file, row.names = FALSE)
}


# Script with commands
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  
output_dir <- args[2]  
spacer <- args[3]     


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

process_file(input_file, output_dir, spacer)