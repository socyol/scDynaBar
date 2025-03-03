##############################################
###     Gastruloids ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

#_______________________________________
# 2. Load metadata file Gastruloids
# ======================================
input_mt <- "data/metadatas/"
metadata <- read.csv(file.path(input_mt, "/metadata_gastruloids.csv")) 
metadata <- na.omit(metadata)


# Selection of cell populations with more than 50 cells
selected_celltypes <- c("Caudal epiblast", "Caudal Mesoderm", "Intermediate mesoderm", "NMP", "Paraxial mesoderm", "Somitic mesoderm")
filtered_metadata <- metadata %>%
  filter(celltype %in% selected_celltypes)
filtered_metadata<- filtered_metadata[filtered_metadata$Coverage >= 4 & filtered_metadata$PassQC == "TRUE" & filtered_metadata$GFP_Status == "Positivo", ]


#_______________________________________
# 3.  Plots
# ======================================

#   Figure 5c :  Boxplot % original sequence in each celltype
# ================================================================

boxplot_metadata_feature <- function(data, feature_name) {
  plot_data <- data.frame(
    Feature = data[[feature_name]]/100,
    CellType = factor(data$celltype)  
  )
  
  # Kruskal-Wallis
  kruskal_test <- kruskal.test(Feature ~ CellType, data = plot_data)
  p_value <- kruskal_test$p.value
  
  purple_colors <- c("#DADAEB","#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#4A1486")
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
  
  # return(kruskal_test)
  return(gg)
}

boxplot_metadata_feature(filtered_metadata, "Uncuts")



#   Figure 5d :   Boxplot showing the distribution of the proportion of original barcodes across different time points
# ================================================================
timecourse <- read.csv(file.path(input_mt, "/metadata_timecourse.csv"))

filtered_gast <- metadata[metadata$Coverage >= 4, ]

# split by day
metadata_day0 <- subset(timecourse, Day == 0 & Coverage >= 4)
metadata_day4 <- subset(timecourse, Day == 4  & Coverage >=  4)
metadata_day10 <- subset(timecourse, Day == 10 & Coverage >= 4)

uncuts_day6 <- sort(filtered_gast$Uncuts)
uncuts_day0 <- metadata_day0$uncut 
uncuts_day4 <- metadata_day4$uncut 
uncuts_day10 <- metadata_day10$uncut 

# day 0:
mean_day0 <- mean(uncuts_day0, na.rm = TRUE)
sd_day0 <- sd(uncuts_day0, na.rm = TRUE)
lower_limit <- mean_day0 - 2 * sd_day0
upper_limit <- mean_day0 + 2 * sd_day0
uncuts_day0_filtered <- uncuts_day0[uncuts_day0 >= lower_limit & uncuts_day0 <= upper_limit]

remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)  
  Q3 <- quantile(x, 0.75, na.rm = TRUE)  
  IQR <- Q3 - Q1                         
  lower_bound <- Q1 - 1.5 * IQR          
  upper_bound <- Q3 + 1.5 * IQR          
  x[x >= lower_bound & x <= upper_bound] 
}

uncuts_day0_filtered <- remove_outliers(uncuts_day0_filtered)  # Day 0 ya filt
#       --- 

colors <- c("Day 0" = "#FDDC6C", "Day 4" = "lightsalmon3", "Day 6" = "#DECBE4", "Day 10" = "#79CDCD")


data_boxplot <- data.frame(
  uncuts = c(uncuts_day0_filtered, uncuts_day4, uncuts_day6, uncuts_day10),
  Day = factor(rep(c("Day 0", "Day 4", "Day 6", "Day 10"),
                   times = c(length(uncuts_day0_filtered), length(uncuts_day4), 
                             length(uncuts_day6), length(uncuts_day10))),
               levels = c("Day 0", "Day 4", "Day 6", "Day 10")) 
)

ggplot(data = data_boxplot, aes(x = Day, y = uncuts, fill = Day)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Day)) +  
  geom_jitter(aes(color = Day), width = 0.2, alpha = 1) +     
  scale_fill_manual(values = colors) +  
  scale_color_manual(values = colors) + 
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +  
  labs(title = "",                      
       x = "Day",                        
       y = "uncuts") +               
  theme_minimal()    



#   Supplementary Figure 10A: Boxplot showing the barcode diversity score for each cell grouped by cell type
# ================================================================
boxplot_metadata_feature <- function(data, feature_name) {
  if (!feature_name %in% colnames(data)) {
    stop("El dataframe no contiene la columna especificada.")
  }
  
  plot_data <- data.frame(
    Feature = data[[feature_name]],
    CellType = factor(data$celltype)  
  )
  
  kruskal_test <- kruskal.test(Feature ~ CellType, data = plot_data)
  p_value <- kruskal_test$p.value

  purple_colors <- c("#DADAEB","#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#4A1486")
  gg <- ggplot(plot_data, aes(x = CellType, y = Feature, fill = CellType)) +
    geom_boxplot(alpha = 0.5, width = 0.5, aes(color = CellType), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.25, aes(color = CellType), size = 1) +
    xlab("") +
    ylab("Diversity") + 
    theme_minimal() +
    geom_jitter(aes(color = CellType), width = 0.2, alpha = 1) +
    scale_fill_manual(values = purple_colors) +
    scale_color_manual(values = purple_colors) +
    # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 5, color = "white", fill = "black") +
    theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"),
          axis.title = element_text(size = 19, color = "black"),
          axis.text = element_text(size = 19, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(color = "gray85", size = 0.5),
          panel.grid.minor = element_line(color = "gray85", size = 0.25),
          legend.position = "none") +
    coord_cartesian(ylim = c(0, 0.45))
  
  # return(kruskal_test)
  return(gg)
}

boxplot_metadata_feature(filtered_metadata, "Diversity")



#   Supplementary Figure 10B: Boxplot showing the barcode diversity score across different time points in single cells 
# ================================================================
diversity_day6 <- sort(filtered_metadata$Diversity) 
diversity_day0 <- metadata_day0$Diversity
diversity_day4 <- metadata_day4$Diversity 
diversity_day10 <- metadata_day10$Diversity

## remove outliers:
diversity_day0 <- remove_outliers(diversity_day0) 
mean_day0 <- mean(diversity_day0, na.rm = TRUE)
sd_day0 <- sd(diversity_day0, na.rm = TRUE)
lower_limit <- mean_day0 - 2 * sd_day0
upper_limit <- mean_day0 + 2 * sd_day0

diversity_day0_filtered <- diversity_day0[diversity_day0 >= lower_limit & diversity_day0 <= upper_limit]


data_boxplot_filtered <- data.frame(
  diversity = c(diversity_day0_filtered, diversity_day4, diversity_day6, diversity_day10),
  Day = factor(rep(c("Day 0", "Day 4", "Day 6", "Day 10"),
                   times = c(length(diversity_day0_filtered), length(diversity_day4), 
                             length(diversity_day6), length(diversity_day10))),
               levels = c("Day 0", "Day 4", "Day 6", "Day 10"))  
)


kruskal_test <- kruskal.test(diversity ~ Day, data = data_boxplot_filtered)

ggplot(data = data_boxplot_filtered, aes(x = Day, y = diversity, fill = Day)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Day)) +  
  geom_jitter(aes(color = Day), width = 0.2, alpha = 1) +     
  scale_fill_manual(values = colors) +  
  scale_color_manual(values = colors) + 
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +  # Marcar la media
  labs(title = "",                      
       x = "Day",                        
       y = "Diversity") +               
  theme_minimal()+
  coord_cartesian(ylim = c(0, 0.45))+
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),   
    axis.text.y = element_text(size = 14, face = "bold"),    
    axis.title.x = element_text(size = 16, face = "bold"),   
    axis.title.y = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14),                   
    legend.title = element_text(size = 16, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  )