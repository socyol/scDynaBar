##############################################
###      Zscan4 ANALYSIS >>> Plots     ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

#_______________________________________
# 2. Load metadata file 
# ======================================
input_mt <- "data/metadatas/"
metadata <- read.csv(file.path(input_mt, "metadata_zscan4.csv")) 
metadata <- na.omit(metadata)

input_so <- "data/seurat_objects/"
sce_zscan4 <- readRDS(file.path(input_so, "so_zscan4_r.rds"))


#_______________________________________
# 3.  Compute expression zscan4 genes
# ======================================

# Compute the expression of zscan4 in cells without cluster 1
sce_cl02 <- subset(sce_zscan4, TSNE_cluster != 1)
all.genes <- rownames(sce_cl02) 

genes_zscan <- c("Zscan4d", "Zscan4c", "Zscan4f", "Zscan4b", "Zscan4e")
for (gene_name in genes_zscan) {
  gene_ids <- grep(paste0("^", gene_name), all.genes, value = TRUE)
  gene_expression <- GetAssayData(object = sce_cl02, layer = "data")[gene_ids, ]
  sce_cl02@meta.data[[tolower(gene_name)]] <- gene_expression
}

# Merge metadata with the seurat object
metadata <- subset(metadata, Coverage >= 4)
colnames_to_keep <- intersect(colnames(sce_cl02), metadata$CellID)
sce_zscan4_reads <- sce_cl02[, colnames_to_keep]

sce_zscan4_reads@meta.data$CellID <- rownames(sce_zscan4_reads@meta.data)
sce_zscan4_reads@meta.data <- merge(sce_zscan4_reads@meta.data, metadata, by = "CellID")
rownames(sce_zscan4_reads@meta.data) <- sce_zscan4_reads@meta.data$CellID


#_______________________________________
# 4.  Plots
# ======================================

#   Figure 4c :  HEATMAP Plot
# ===============================

# Define totipotency  genes
functional_genes <- c("Zscan4d", "Usp17lc", "Usp17la", "Usp17le", "Usp17lb", "Usp17ld", "Tcstv3", "Zscan4c", "Zscan4e", "Zscan4b", "Zscan4f", "Sp110", "Ankrd22", "Abcb5", "Cacna1s", "Popdc3", "Dkk1", "Ccl3", "Lgals4", "Pramel7", "Usp9y", "Zfp352", "Adh7", "Phf11a", "Slfn2")

genes_heatmap <- functional_genes 
sce_heatmap <- subset(sce_zscan4_reads, TSNE_cluster == 0 | TSNE_cluster == 2)

# If needed, reduce the number of cells in cluster 0
cells_cluster_0 <- subset(sce_zscan4_reads, TSNE_cluster == 0)
set.seed(123)
random_cells_cluster_0 <- cells_cluster_0@meta.data[sample(rownames(cells_cluster_0@meta.data), 200), ]
cells_cluster_2 <- subset(sce_zscan4_reads, TSNE_cluster == 2)

sce_heatmap <- subset(sce_zscan4_reads, cells = c(rownames(random_cells_cluster_0), rownames(cells_cluster_2@meta.data)))

# Generate heatmap matrix
sce_zscan4_heatmap <- sce_heatmap[genes_heatmap,]
expression_matrix <- GetAssayData(sce_zscan4_heatmap, layer = "data")
dense_expression_matrix_toti <- as.matrix(expression_matrix)

cluster_info <- sce_zscan4_heatmap@meta.data$TSNE_cluster
annotation_col<- data.frame(TSNE_cluster = cluster_info)
rownames(annotation_col) <- colnames(dense_expression_matrix_toti)  

# Define colors
cluster_colors <- list(TSNE_cluster = c("0" = "#7AC5CD", "2" = "#FFA07A")) 
heatmap_colors <- colorRampPalette( c( "#467CBD", "#AAC5E5", "#FBE8C5", "#F2DAAC" , "#F07837", "#EA4A2E"))(20)
col_annotation <- HeatmapAnnotation(df = annotation_col, col = cluster_colors)

ht1 <- Heatmap(dense_expression_matrix_toti,
               name = "Expression",
               col = heatmap_colors, 
               show_row_names = TRUE, 
               show_column_names = FALSE,
               cluster_rows = FALSE,  
               cluster_columns = TRUE, 
               clustering_distance_rows = "euclidean",
               # clustering_method_rows = "complete",
               top_annotation = col_annotation,
               row_title  = "Totipotent related genes")  
ht1


#    Figure 4d: Zscan4c vs Diversity colored by Uncuts
# ====================================================
metadatap <- sce_zscan4_reads@meta.data
metadatap <- subset(metadatap, Coverage >= 8)

palete <- brewer.pal(n = 9, "RdPu")[3:9]
ggplot(metadatap, aes(x = Diversity, y = zscan4c, color = Uncuts, fill = Uncuts)) +
  geom_point(size = 4, shape = 16) + 
  scale_color_gradientn(colours = palete, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete, name = "% Uncuts") +
  labs(x = "Diversity", y = "Log Zscan4c expression") +
  theme_classic()+
  ggtitle("Zscan4+ mESC")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 19, color = "black", face = "bold"), 
    axis.title = element_text(size = 19, color = "black"), 
    axis.text = element_text(size = 19, color = "black"))




#   Supplementary Figure 8: Boxplot showing the barcoding cassette coverage in 2C-like cells and pluripotent cells
# ====================================================

sce_zscan4_reads@meta.data[["TSNE_cluster"]] <- factor(sce_zscan4_reads@meta.data[["TSNE_cluster"]], levels = c("1", "2", "0"))# c("1", "2", "0")
levels(sce_zscan4_reads@meta.data$TSNE_cluster) <- c("Control", "Totipotents\nZscan4+", "Pluripotents\nZscan4-")

zscan4_meta_1 <- subset(sce_zscan4_reads@meta.data, TSNE_cluster != "Control")
zscan4_meta_1$TSNE_cluster <- factor(
  zscan4_meta_1$TSNE_cluster,
  levels = c("Totipotents\nZscan4+", "Pluripotents\nZscan4-"))

levels(zscan4_meta_1$TSNE_cluster)
table(zscan4_meta_1$TSNE_cluster)

accent_colors <- brewer.pal(8, "Accent")
colors <- c(accent_colors[3], accent_colors[1])

ggplot(zscan4_meta_1, aes(x = TSNE_cluster, y = Coverage, fill = TSNE_cluster)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = TSNE_cluster), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.25, aes(color = TSNE_cluster), size = 1) + # Más grueso
  xlab("") +
  ylab("Coverage gRNA locus") + # "%Original sequence" "Diversity per cell" Mean length (bp)
  theme_minimal() +
  geom_jitter(aes(color = TSNE_cluster), width = 0.2, alpha = 1) +
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "black", fill = "black") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ylim(0,300)+
  theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"), # Título más oscuro y grande
        axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
        axis.text = element_text(size = 19, color = "black"),
        panel.grid.major = element_line(color = "gray85", size = 0.5),  # Adjusts major grid lines
        panel.grid.minor = element_line(color = "gray85", size = 0.25), # Texto de ejes más oscuro
        legend.position = "none")



## kruskal_test
kruskal.test(Coverage ~ TSNE_cluster, data = zscan4_meta_1)

## wilcoxon test
group1 <- subset(zscan4_meta_1, TSNE_cluster == "Pluripotents\nZscan4-")
group2 <- subset(zscan4_meta_1, TSNE_cluster == "Totipotents\nZscan4+")

wilcox_test <- wilcox.test(group1$Coverage, group2$Coverage)
p_value <- wilcox_test$p.value
print(wilcox_test$p.value)


#   Supplementary Figure 12b: QC Illumina Zscan4 exp
# ====================================================

#        🚨   To do this you have to previously run script 2_sc_Barcode-AlelleSequences_analysis.R  🚨
#             --------------------------------------------------------------------------------------

## 1. Load the reads table result
input_results <- "data/results" # in case you've saved the reads_Table data frames
reads_table_zscan4 <- read.csv(file.path(input_results, "NAME_METADATA")) 

## 2. Select the control group of the filtered seurat object
cells_cl_1 <- rownames(sce_zscan4@meta.data)[sce_zscan4@meta.data$TSNE_cluster == 1] 
filtered_barcodes <- reads_table_zscan4[reads_table_zscan4$Cell_ID %in% cells_cl_1, ]


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
  labs(title = "Diversity distribution Control Cluster (1)",
       x = "Illumina Quality Scores ranges",
       y = "Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

