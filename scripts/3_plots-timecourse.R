##############################################
###      TIMECOURSE ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("settings.R")

#_______________________________________
# 2. Load metadata file 
# ======================================

input_mt <- "data/metadatas/"
metadata <- read.csv(file.path(input_mt, "metadata_timecourse.csv")) 
metadata <- na.omit(metadata)

#  Figure 3c: Boxplot % Original barcode sequence
# ====================================================
# Define the colours
colors <- c("0" = "darkgoldenrod1", "4" = "red", "10" = "#00CDCD")

ggplot(data = metadata, aes(x = factor(Day), y = Uncuts, fill = factor(Day))) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = factor(Day)), outlier.shape = NA) +  
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "",
       x = "Day",
       y = "% Original sequence") +
  theme_minimal()


#  Figure 3e: Boxplot mean sequence diversity
# ====================================================

ggplot(data = metadata, aes(x = factor(Day), y = Diversity, fill = factor(Day))) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = factor(Day)), outlier.shape = NA) + 
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "",
       x = "Day",
       y = "Diversity per cell") +
  theme_minimal()

#  Figure 3f: Boxplotmean length
# ====================================================

ggplot(data = metadata, aes(x = factor(Day), y = Mean_length, fill = factor(Day))) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = factor(Day)), outlier.shape = NA) + 
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(title = "",
       x = "Day",
       y = "Mean length (bp)") +
  theme_minimal()



#   Supplementary Figure 5a: QC Illumina Timecourse
# ====================================================

## 1. Load the reads table result
input_results <- "data/results" # in case you've saved the reads_Table data frames
reads_table_timecourse <- read.csv(file.path(input_results, "NAME_METADATA")) 
View(datostimecourse)


## 2. Select cells IDs day 0
cell_ids_dia_0 <- metadata$CellID[metadata$Day == 0]
cell_ids_dia_0<-na.omit(cell_ids_dia_0)
filtered_barcodes<- reads_table_timecourse[reads_table_timecourse$Cell_ID %in% cell_ids_dia_0, ]

# modified code for processing the data, without filtering
process_data <- function(data, output_path = "") {
  
  # 1) Remove rows with Length < 4 or > 35 and remove NA values
  data <- subset(data, Length > 4 & Length < 35)
  data <- na.omit(data)
  
  # 2) Add a "Diversity" column based on the negative of "Alignment_score"
  data$Diversity <- -data$Alignment_score
  data$Diversity <- (data$Diversity - min(data$Diversity)) / 
    (max(data$Diversity) - min(data$Diversity))
  
  return(data)
}

data <- process_data(filtered_barcodes)

## 3. Add a column with the ranges
data$Quality_Scores_Range <- cut(data$Quality_Scores, 
                                 breaks = seq(16, 38, by = 2), 
                                 include.lowest = TRUE, 
                                 right = FALSE)

## 4. plot of "beans" or boxplots
cell_counts <- data %>% ## ---------> to add number of cells
  group_by(Quality_Scores_Range) %>%
  summarise(count = n())

ggplot(data = data, aes(x = Quality_Scores_Range, y = Diversity)) +
  geom_boxplot() +
  geom_text(data = cell_counts, aes(x = Quality_Scores_Range, y = max(data$Diversity) + 0.1, label = count), 
            vjust = -0.5, size = 3) + 
  labs(title = "Diversity distribution Day 0",
       x = "Illumina Quality Scores ranges",
       y = "Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

