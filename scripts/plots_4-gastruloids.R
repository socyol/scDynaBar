##############################################
###     Gastruloids ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

#_______________________________________
# 2. Load metadata file 
# ======================================
input_mt <- "data/metadatas/"
input_mt <- "Documents/TFM/paper_reduction/data/metadatas"

metadata <- read.csv(file.path(input_mt, "/metadata_gastruloids.csv")) 
metadata <- na.omit(metadata)

#_______________________________________
# 3.  Plots
# ======================================

#   Figure 5c :  Boxplot % original sequence in each celltype
# ================================================================
# selecion of cell populaions with more than 50 cells
selected_celltypes <- c("Caudal epiblast", "Caudal Mesoderm", "Intermediate mesoderm", "NMP", "Paraxial mesoderm", "Somitic mesoderm")
filtered_metadata <- metadata %>%
  filter(celltype %in% selected_celltypes)
# 
# 
# Q1 <- quantile(filtered_metadata$Uncuts, 0.25, na.rm = TRUE)
# Q3 <- quantile(filtered_metadata$Uncuts, 0.75, na.rm = TRUE)
# IQR <- Q3 - Q1
# 
# lower_bound <- Q1 - 1.5 * IQR
# upper_bound <- Q3 + 1.5 * IQR
# 
# 
# filtered_metadata <- filtered_metadata %>%
#   filter(Uncuts >= lower_bound & Uncuts <= upper_bound)
# 
# ggplot(filtered_metadata, aes(x = factor(celltype), y = Uncuts, fill = factor(celltype))) +
#   geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) + 
#   geom_jitter(color = "black", size = 1.5, alpha = 0.5, width = 0.2) +  
#   scale_fill_brewer(palette = "Purples") + 
#   labs(title = "", x = "Cell Type", y = "% Original sequence") +
#   theme_minimal() +
#   theme(legend.position = "none")  
# 


boxplot_metadata_feature <- function(data, feature_name) {
  if (!feature_name %in% colnames(data)) {
    stop("El dataframe no contiene la columna especificada.")
  }
  
  plot_data <- data.frame(
    Feature = data[[feature_name]],
    CellType = factor(data$celltype)  
  )
  
  for (celltype in unique(plot_data$CellType)) {
    subset_data <- plot_data[plot_data$CellType == celltype, ]
    Q1 <- quantile(subset_data$Feature, 0.25, na.rm = TRUE)
    Q3 <- quantile(subset_data$Feature, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    plot_data <- plot_data[!(plot_data$CellType == celltype & 
                               (plot_data$Feature < lower_bound | plot_data$Feature > upper_bound)), ]
  }

  purple_colors <- brewer.pal(n = length(unique(plot_data$CellType)), name = "Purples")
  gg <- ggplot(plot_data, aes(x = CellType, y = Feature, fill = CellType)) +
    geom_boxplot(alpha = 0.5, width = 0.5, aes(color = CellType), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.25, aes(color = CellType), size = 1) +
    xlab("") +
    ylab("%Original sequence") + 
    theme_minimal() +
    geom_jitter(aes(color = CellType), width = 0.2, alpha = 1) +
    scale_fill_manual(values = purple_colors) +
    scale_color_manual(values = purple_colors) +
    theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"),
          axis.title = element_text(size = 19, color = "black"),
          axis.text = element_text(size = 19, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(color = "gray85", size = 0.5),
          panel.grid.minor = element_line(color = "gray85", size = 0.25),
          legend.position = "none")
  
  return(gg)
}

boxplot_metadata_feature(filtered_metadata, "Uncuts")


#   Figure 5d :   Boxplot showing the distribution of the proportion of original barcodes across different time points
# ================================================================
timecourse <- read.csv(file.path(input_mt, "/metadata_timecourse.csv"))
filtered_gast <- metadata[metadata$Coverage > 4, ]

# split by day
metadata_day0 <- subset(timecourse, Day == 0 & Coverage >4)
metadata_day4 <- subset(timecourse, Day == 4  & Coverage > 4)
metadata_day10 <- subset(timecourse, Day == 10 & Coverage >4)

uncuts_day6 <- sort(filtered_gast$Uncuts)
uncuts_day0 <- metadata_day0$uncut *100
uncuts_day4 <- metadata_day4$uncut * 100
uncuts_day10 <- metadata_day10$uncut * 100

colors <- c("Day 0" = "#FDDC6C", "Day 4" = "lightsalmon3", "Day 6" = "#DECBE4", "Day 10" = "#79CDCD")

data_boxplot <- data.frame(
  diversity = c(diversity_day0, diversity_day4, diversity_day6, diversity_day10),
  Day = factor(rep(c("Day 0", "Day 4", "Day 6", "Day 10"),
                   times = c(length(diversity_day0), length(diversity_day4), 
                             length(diversity_day6), length(diversity_day10))),
               levels = c("Day 0", "Day 4", "Day 6", "Day 10"))  # Establecer el orden de los niveles
)

data_boxplot <- data.frame(
  uncuts = c(uncuts_day0, uncuts_day4, uncuts_day6, uncuts_day10),
  Day = factor(rep(c("Day 0", "Day 4", "Day 6", "Day 10"),
                   times = c(length(uncuts_day0), length(uncuts_day4), 
                             length(uncuts_day6), length(uncuts_day10))),
               levels = c("Day 0", "Day 4", "Day 6", "Day 10")) 
)

ggplot(data = data_boxplot, aes(x = Day, y = uncuts, fill = Day)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Day)) +  
  geom_jitter(aes(color = Day), width = 0.2, alpha = 1) +     
  scale_fill_manual(values = colors) +  
  scale_color_manual(values = colors) + 
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +  # Marcar la media
  labs(title = "",                      
       x = "Day",                        
       y = "uncuts") +               
  theme_minimal()    


