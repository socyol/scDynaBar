##############################################
###      TIMECOURSE ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

#_______________________________________
# 2. Load metadata file 
# ======================================
metadata <- read.csv("/home/ibmb-ihh-02/Documents/TFM/scDynaBar/data/metadatas/metadata_timecourse.csv") # IRENE METADATA!

input_mt <- "data/metadatas/"
metadata <- read.csv(file.path(input_mt, "metadata_timecourse.csv")) 
metadata <- na.omit(metadata)

#  Figure 3c: Boxplot % Original barcode sequence
# ====================================================
# Define the colours
if (is.factor(metadata$uncut)) {
  metadata$uncut <- as.numeric(as.character(metadata$uncut))
} else if (is.character(metadata$uncut)) {
  metadata$uncut <- as.numeric(metadata$uncut)
}

metadata$Day <- factor(metadata$Day, levels = c(0, 4, 10))

ggplot(data = metadata, aes(x = factor(Day), y = uncut, fill = factor(Day))) +
  geom_boxplot(width = 0.7, alpha = 0.5, color = c("#fed976",  "#F8766D", "#00BFC4" ), outlier.shape = NA) +
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = c("#fed976",  "#F8766D", "#00BFC4" ) ) +  # Custom fill colors
  scale_color_manual(values = c("#fed976",  "#F8766D", "#00BFC4" )) +  # Match jitter color to boxplot fill
  theme_minimal() +
  labs(title = "",
       x = "Day",
       y = "% Original sequence") 

#  Figure 3b: Boxplotmean intact PAM
# ====================================================
ggplot(data = metadata, aes(x = factor(Day), y = pam, fill = factor(Day))) +
  geom_boxplot(width = 0.7, alpha = 0.5, color = c("#fed976",  "#F8766D", "#00BFC4" ), outlier.shape = NA) +
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = c("#fed976",  "#F8766D", "#00BFC4" ) ) +  # Custom fill colors
  scale_color_manual(values = c("#fed976",  "#F8766D", "#00BFC4" )) +  # Match jitter color to boxplot fill
  theme_minimal() +
  labs(title = "",
       x = "Day",
       y = "% Original sequence") 

#  Figure 3e: Boxplot mean sequence diversity
# ====================================================
ggplot(data = metadata, aes(x = factor(Day), y = Diversity, fill = factor(Day))) +
  geom_boxplot(width = 0.7, alpha = 0.5, color = c("#fed976",  "#F8766D", "#00BFC4" ), outlier.shape = NA) +
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = c("#fed976",  "#F8766D", "#00BFC4" ) ) +  # Custom fill colors
  scale_color_manual(values = c("#fed976",  "#F8766D", "#00BFC4" )) +  # Match jitter color to boxplot fill
  theme_minimal() +
  labs(title = "",
       x = "Day",
       y = "Diversity per cell") 


#  Figure 3f: Boxplotmean length
# ====================================================
ggplot(data = metadata, aes(x = factor(Day), y = length, fill = factor(Day))) +
  geom_boxplot(width = 0.7, alpha = 0.5, outlier.shape = NA, aes(color = factor(Day))) +
  geom_jitter(aes(color = factor(Day)), width = 0.2, alpha = 1) +
  scale_fill_manual(values = c("#fed976",  "#F8766D", "#00BFC4")) +  # Custom fill colors
  scale_color_manual(values = c("#fed976",  "#F8766D", "#00BFC4")) +  # Match jitter color to boxplot fill
  theme_minimal() +
  labs(title = "",
       x = "Day",
       y = "Mean length( per cell (bp)")


#   Supplementary Figure 7: Number of unique barcode sequences per cell at Day 0, Day 4 and Day 10
# ====================================================
input_mt <- "data/barcode_sequences/"
data_day0 <- read.csv(file.path(input_mt, "barcode_sequences_tc_day0_filt.csv")) 
data_day4 <- read.csv(file.path(input_mt, "barcode_sequences_tc_day4_filt.csv")) 
data_day10 <- read.csv(file.path(input_mt, "barcode_sequences_tc_day10_filt.csv")) 

data_day0_filtered <- data_day0 %>%
  filter(
    !((Length == 24 & Original == "NO") | 
        (Length == 25 & grepl("G{4}$", Match))) 
  )

unique_barcodes_day0 <- calculate_unique_barcodes_per_cell(data_day0_filtered)
unique_barcodes_day4 <- calculate_unique_barcodes_per_cell(data_day4)
unique_barcodes_day10 <- calculate_unique_barcodes_per_cell(data_day10)
unique_barcodes_day0$Day <- 0
unique_barcodes_day4$Day <- 4
unique_barcodes_day10$Day <- 10

combined_unique_barcodes <- rbind(unique_barcodes_day0, unique_barcodes_day4, unique_barcodes_day10)
ggplot(combined_unique_barcodes, aes(x = factor(Day), y = Unique_Barcodes, fill = factor(Day))) +
  geom_boxplot(color = "#474747", alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(
    values = c("#fed976",  "#F8766D", "#00BFC4" ) 
  )+
  labs(
    x = "Day",
    y = "Number of Unique Barcodes per Cell",
    fill = "Day",
    title = "Distribution of Unique Barcodes per Cell by Day"
  )+
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),    
    axis.text.y = element_text(size = 14, face = "bold"),   
    axis.title.x = element_text(size = 16, face = "bold"),   
    axis.title.y = element_text(size = 16, face = "bold"),   
    legend.text = element_text(size = 14),                   
    legend.title = element_text(size = 16, face = "bold"),   
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  )

#   Supplementary Figure 12a: QC Illumina Timecourse
# ====================================================

#        🚨   To do this you have to previously run script 2_sc_Barcode-AlelleSequences_analysis.R  🚨
#             --------------------------------------------------------------------------------------

## 1. Load the reads table result
input_results <- "data/results" # in case you've saved the reads_Table data frames
reads_table_timecourse <- read.csv(file.path(input_results, "NAME_METADATA")) 


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