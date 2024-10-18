##############################################
###        Bulk analysis Figures          ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("settings.R")

#_______________________________________
# 2. Load metadata file 
# ======================================
input_bs <- "data/metadatas/"
metadata <- read.csv(file.path(input_bs, "metadata_bulk.csv")) 

# select only replicate 1
data_r1 <- metadata %>%
  filter(Replicate == "Replicate1", Coverage > 200)


# statistic tests
jt_test_uncuts <- JonckheereTerpstraTest(data_r1$p.uncut, data_r1$Day)
jt_test_diversity <- JonckheereTerpstraTest(data_r1$Diversity, data_r1$Day)


#_______________________________________
# 3.  Plots
# ======================================

#   Figure 1d :  Proportion of original barcodes over time across different gRNAs
# =======================================================================================
ggplot(data_r1, aes(x = as.factor(Day), y = p.uncut, color = System, group = System)) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) +  
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  scale_shape_manual(values = c(16, 17)) + 
  labs(title = "",
       x = "",
       y = "% Original sequence") +
  facet_wrap(~gRNA) +  
  theme_minimal() +
  theme(legend.position = "right",  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))  

#   Figure 1e :  Proportion of diversity over time across different gRNAs
# =======================================================================================
ggplot(data_r1, aes(x = as.factor(Day), y = Diversity, color = System, group = System)) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) + 
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  scale_shape_manual(values = c(16, 17)) +  
  labs(title = "",
       x = "",
       y = "Barcode Diversity") +
  facet_wrap(~gRNA) +  
  theme_minimal() +
  theme(legend.position = "right",  # Ajustar la posición de la leyenda
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) 


# Figure 1f :  Proportion of original barcodes over time for gRNAs 21 and 26-bp spacers
# =======================================================================================
data_r1_bp <- data_r1 %>%
  mutate(spacer = case_when(
    gRNA %in% c("g1", "g2", "g3", "g4", "g5") ~ "21bp",
    gRNA %in% c("g6", "g9") ~ "26bp",
    TRUE ~ NA_character_  # Para otros gRNAs, puedes poner NA o dejarlo en blanco
  ))
ggplot(data_r1_bp, aes(x = factor(Day), y = p.uncut)) +
  geom_boxplot() +
  # scale_fill_manual(values = c("#2ca25f",  "#e6550d")) + 
  facet_wrap(~ spacer) +  # Crear un gráfico separado por cada condición
  labs(x = "Day", y = "% Original sequence", title = "") +
  theme_minimal()


#   Figure 2d : . Barcode diversity over time
# =======================================================================================
ggplot(data_r1, aes(x = factor(Day), y = Diversity, fill = System)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#2ca25f",  "#e6550d")) + 
  facet_wrap(~ System) +  # Crear un gráfico separado por cada condición
  labs(x = "Day", y = "Barcode Diversity", title = "") +
  theme_minimal()


#   Supplementary Figure 3 : Average of indels per read across time points
# =======================================================================================
ggplot(data_r1, aes(x = System, y = Indels, fill = System)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  labs(title = "Proportion of indels in Cas9 vs BE3",
       x = "System",
       y = "% Indels") +
  theme_minimal()

wilcox_test <- wilcox.test(Indels ~ System, data = data_r1)




