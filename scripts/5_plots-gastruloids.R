##############################################
###     Gastruloids ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("settings.R")

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
selected_celltypes <- c("Caudal epiblast", "Caudal Mesoderm", "Intermediate mesoderm", "NMP", "Paraxial mesoderm", "Somitic mesoderm")
filtered_metadata <- metadata %>%
  filter(celltype %in% selected_celltypes)


Q1 <- quantile(filtered_metadata$Uncuts, 0.25, na.rm = TRUE)
Q3 <- quantile(filtered_metadata$Uncuts, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR


filtered_metadata <- filtered_metadata %>%
  filter(Uncuts >= lower_bound & Uncuts <= upper_bound)

ggplot(filtered_metadata, aes(x = factor(celltype), y = Uncuts, fill = factor(celltype))) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) + 
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, width = 0.2) +  
  scale_fill_brewer(palette = "Purples") + 
  labs(title = "", x = "Cell Type", y = "% Original sequence") +
  theme_minimal() +
  theme(legend.position = "none")  



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

