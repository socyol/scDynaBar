##############################################
###        Bulk analysis Figures          ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

library(clinfun)
library(DescTools)
library(tidyr)
library(gplots)
library(viridisLite)

#_______________________________________
# 2. Load metadata file 
# ======================================
input_bs <- "data/metadatas/"
metadata <- read.csv(file.path(input_bs, "metadata_bulk_df.csv")) 

metadata <- read.csv("/home/ibmb-ihh-02/Documents/TFM/scDynaBar/data/metadatas/metadata_bulk_df.csv")
# select only replicate 1
metadata <- metadata %>%
  filter(X != 130)


# statistic tests
data_r1 <- metadata %>%
  filter(replicate == "Replicate1", Coverage > 200)
jt_test_uncuts <- JonckheereTerpstraTest(data_r1$p.uncut, data_r1$Day)
jt_test_diversity <- JonckheereTerpstraTest(data_r1$Diversity, data_r1$Day)


#_______________________________________
# 3.  Plots
# ======================================

#   Figure 1c :  . Proportion of barcodes based on their mutational profile: active, inactive); and original 
# =======================================================================================
df<-metadata
df$p.active.uncut.PAM <- df$active.UncutPAM.qualiy/df$Coverage
df$p.inactive.cut  <- df$inactive.PAM.quality/df$Coverage
df <-  df[df$Coverage > 200,]

# aggregate(df$p.uncut.quality[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.uncut
# aggregate(df$p.active.uncut.PAM[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.active.uncutPAM
# aggregate(df$p.inactive.cut[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.inactive
# df2 <- data.frame( Day= rep(mean.day.inactive[,1], 3), value= c( mean.day.inactive[,2],mean.day.active.uncutPAM[,2],mean.day.uncut[,2] ),
#                    activity= c( rep("inactive", 6), rep( "activeUncutPAM", 6), rep( "uncut", 6)))
# act <- factor(df2$activity, levels=c(  "uncut", "activeUncutPAM", "inactive" ))
# ggplot(data=df2, aes(x=as.numeric(as.vector(Day)), y=value, fill= act)) +
#   geom_area() +   scale_fill_brewer(palette="Blues") + theme_classic()


# -----------------
#       Cas 9
# -----------------
df_Cas9 <- df[df$System == "Cas9", ]

mean.day.uncut_Cas9 <- aggregate(df_Cas9$p.uncut[df_Cas9$replicate == "Replicate1"], 
                                 list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate1"]), FUN = mean)
mean.day.active_Cas9 <- aggregate(df_Cas9$p.active[df_Cas9$replicate == "Replicate1"], 
                                  list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate1"]), FUN = mean)
mean.day.inactive_Cas9 <- aggregate(df_Cas9$p.inactive[df_Cas9$replicate == "Replicate1"], 
                                    list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate1"]), FUN = mean)


df2_Cas9 <- data.frame(
  Day = rep(mean.day.uncut_Cas9[, 1], 3),
  value = c(mean.day.uncut_Cas9[, 2], mean.day.active_Cas9[, 2], mean.day.inactive_Cas9[, 2]),
  activity = c(rep("Original", nrow(mean.day.uncut_Cas9)), 
               rep("Active", nrow(mean.day.active_Cas9)), 
               rep("Inactive", nrow(mean.day.inactive_Cas9)))
)

df2_Cas9$activity <- factor(df2_Cas9$activity, levels = c("Original", "Active", "Inactive"))

ggplot(data = df2_Cas9, aes(x = as.numeric(as.vector(Day)), y = value, fill = activity)) +
  geom_area() +
  scale_fill_brewer(palette="Oranges")  +
  theme_classic() +
  labs(x = "Day", y = "Proportion of barcodes", fill = "Activity", title = "Cas9 ")+
  theme(plot.title = element_text(hjust = 0.5))

# ------------------
#     BE3
# ------------------
df_BE3 <- df[df$System == "BE3", ]

mean.day.uncut_BE3 <- aggregate(df_BE3$p.uncut[df_BE3$replicate == "Replicate1"], 
                                list(df_BE3$Day.c[df_BE3$replicate == "Replicate1"]), FUN = mean)
mean.day.active_BE3 <- aggregate(df_BE3$p.active[df_BE3$replicate == "Replicate1"], 
                                 list(df_BE3$Day.c[df_BE3$replicate == "Replicate1"]), FUN = mean)
mean.day.inactive_BE3 <- aggregate(df_BE3$p.inactive[df_BE3$replicate == "Replicate1"], 
                                   list(df_BE3$Day.c[df_BE3$replicate == "Replicate1"]), FUN = mean)


df2_BE3 <- data.frame(
  Day = rep(mean.day.uncut_BE3[, 1], 3),
  value = c(mean.day.uncut_BE3[, 2], mean.day.active_BE3[, 2], mean.day.inactive_BE3[, 2]),
  activity = c(rep("Original", nrow(mean.day.uncut_BE3)), 
               rep("Active", nrow(mean.day.active_BE3)), 
               rep("Inactive", nrow(mean.day.inactive_BE3)))
)

df2_BE3$activity <- factor(df2_BE3$activity, levels = c("Original", "Active", "Inactive"))

ggplot(data = df2_BE3, aes(x = as.numeric(as.vector(Day)), y = value, fill = activity)) +
  geom_area() +
  scale_fill_brewer(palette="BuGn") +
  theme_classic() +
  labs(x = "Day", y = "Proportion of barcodes", fill = "Activity", title = "BE3 ")+
  theme(plot.title = element_text(hjust = 0.5))




#   Figure 1d :  Proportion of original barcodes over time across different gRNAs
# =======================================================================================
data_r1 <- df %>%
  filter(replicate == "Replicate1", Coverage > 200)

ggplot(data_r1, aes(x = as.factor(Day.c), y = p.uncut, color = System, group = System)) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) +  
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  scale_shape_manual(values = c(16, 16)) + 
  labs(title = "",
       x = "",
       y = "% Original sequence") +
  facet_wrap(~gRNA) +  
  theme_classic() +
  theme(legend.position = "right",  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))  


#   Figure 1e :  Proportion of diversity over time across different gRNAs
# =======================================================================================
ggplot(data_r1, aes(x = as.factor(Day.c), y = Diversity, color = System, group = System)) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) + 
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  scale_shape_manual(values = c(16, 16)) +  
  labs(title = "",
       x = "",
       y = "Barcode Diversity") +
  facet_wrap(~gRNA) +  
  theme_classic() +
  theme(legend.position = "right",  
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
    TRUE ~ NA_character_  
  ))
ggplot(data_r1_bp, aes(x = factor(Day), y = p.uncut)) +
  geom_boxplot() +
  facet_wrap(~ spacer) +  
  scale_fill_brewer(palette = "Purples",  direction = 1) +
  theme_classic()+
  labs(x = "Day", y = "% Original sequence", title = "")



#   Figure 1g : Percentage of the original nucleotide over time, aligning the spacer relative to the PAM sequence (Cas9 system)
# =======================================================================================

df <-  metadata[metadata$Coverage > 200,]
my_list <- df$Mutationperbase
# Function to split string and convert to numeric
split_and_convert <- function(num) {
  as.numeric(strsplit(as.character(num), "_")[[1]])
}

# Convert the list to a matrix
matrix_result <- do.call(rbind, lapply(my_list, split_and_convert))

# Display the resulting matrix
print(matrix_result)
mutperbase.mat <- do.call(rbind, mutperbase)  
tmp <- matrix_result
rownames(tmp) <- df$Day.c
tmp <-tmp[order(df$Day),]

tmp2 <- tmp[df$System =="Cas9" ,]
tmp3 <- rbind( colMeans(tmp2[ rownames(tmp2) == 0, ]), colMeans(tmp2[ rownames(tmp2) == 2, ] ), colMeans(tmp2[ rownames(tmp2) == 4, ]),
               colMeans(tmp2[ rownames(tmp2) == 8, ]), colMeans(tmp2[ rownames(tmp2) == 10, ]), colMeans(tmp2[ rownames(tmp2) == 18, ]),
               colMeans(tmp2[ rownames(tmp2) == 25, ]), colMeans(tmp2[ rownames(tmp2) == 32, ]))


rampCol2 <- colorRampPalette(c( "red", "#f03b20", "#feb24c" , "#ffffcc","white"))(n = 1000)


heatmap.2( tmp3 , na.color ="#f1a34", Colv = FALSE, Rowv =FALSE ,trace="none" , col=viridis )
color_palette <- viridis(1000)

heatmap.2(tmp3,
          na.color = "#f1a34",
          Colv = FALSE, Rowv = FALSE,
          trace = "none",
          col = color_palette,
          main = title,
          xlab = "Position relative to PAM",
          ylab = "Day",
          margins = c(5, 5))


#   Figure 1h : Percentage of C>T mutations over time
# =======================================================================================
ggplot(df[ df$replicate == "Replicate1",], aes(x=Day.c, y=C.T,  fill=System)) + geom_boxplot( alpha=0.9) +
  theme_classic()  +    ylim(0, 0.1) +
  scale_fill_manual(values=c("#2ca25f",  "#e6550d"))

df_filtered <- df %>%
  filter(replicate == "Replicate1")

# ✅ Crear el gráfico
ggplot(df_filtered, aes(x = as.factor(Day), y = C.T, fill = System)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1, width = 0.6) +
  scale_fill_manual(values = c("Cas9" = "#e74c3c", "BE3" = "#2ecc71")) +  # Colores similares a la figura
  theme_classic() +
  labs(x = "Days", y = "C > T mutations", fill = "System") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )




#   Figure 2b :  Pearson’s correlation coefficient for reproducibility between replicates from days 0 to 18.
# =======================================================================================
## correlation plot
replicate1 <- c("8227", "8238")
replicate2<- c( "8467", "8468" )

df <- metadata %>%
  filter(Coverage > 200)
df1 <- df[df$replicate == "Replicate1" &  df$Coverage > 200,]
df2 <- df[df$replicate == "Replicate2" & df$Coverage > 200,]

data.cor <-c()
data.mean  <-c()
for ( i in  c(1:nrow(df1))) {
  day= df1$Day.c[i]
  gRNA= df1$gRNA[i]
  exp= df1$System[i]
  value <- which(df2$Day.c == day & df2$gRNA== gRNA & df2$System == exp)  
  if (length(value) ==  1)
  {
    V <- c(day, gRNA, exp, df1$p.uncut.quality[i], df2$p.uncut.quality[value])
    
    data.cor <- rbind( data.cor, V)
  }
}

####
####    For the reviewers they ask us to say how many replicates and combinations we have used for the analysis
####    of the Fig 2B. We created a table (excell type) with the conditions, replicates, gRNAs and system
#     . ........................................
data_table <- data.frame(
  Day = as.numeric(data.cor[, 1]),
  gRNA = data.cor[, 2],
  System = data.cor[, 3],
  p_uncuts_Replicate1 = as.numeric(data.cor[, 4]),
  p_uncuts_Replicate2 = as.numeric(data.cor[, 5])
)
num_combinations <- nrow(data_table)
num_combinations <- n_distinct(data_table[, c("Day", "gRNA", "System")]) # 33

num_days <- length(unique(data_table$Day)) # 6
num_replicates <- 2  

cat(" Total unique combinations (gRNA-System-Day):", num_combinations, "\n",
    "Total days analyzed:", num_days, "\n",
    "Total replicates per combination:", num_replicates, "\n")


# ..............................................
replicateA <- as.numeric(data.cor[,4])
replicateB <- as.numeric(data.cor[,5])

data <- data.frame(replicateA, replicateB)
correlation <- cor(replicateA, replicateB)
model <- lm(replicateB ~ replicateA, data = data)
x_seq <- seq(min(replicateA), max(replicateA), length.out = 100)
pred <- predict(model, newdata = data.frame(replicateA = x_seq), interval = "confidence")
shading_data <- data.frame(x = c(x_seq, rev(x_seq)), y = c(pred[, "lwr"], rev(pred[, "upr"])))

# Create the scatter plot with shaded area around the correlation line
ggplot(data, aes(x = replicateA, y = replicateB)) +
  geom_point( size=5, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_polygon(data = shading_data, aes(x = x, y = y), fill = "red", alpha = 0.1) +
  labs( x = "Replicate 1",
        y = "Replicate 2") +
  annotate("text", x = min(data$replicateA), y = max(data$replicateB),
           label = paste("Correlation =", round(correlation, 2)), hjust = 0, vjust = 1) +
  theme_classic()


#   Figure 2c :  . Proportion of barcodes based on their mutational profile: Cas9 and BE3
# =======================================================================================
df <- metadata %>%
  filter(Coverage > 200)

#     Cas9:
# -------------
df_Cas9 <- df[df$System == "Cas9", ]
mean.day.uncut_Cas9 <- aggregate(df_Cas9$p.uncut, list(df_Cas9$Day.c), FUN = mean)
mean.day.active_Cas9 <- aggregate(df_Cas9$p.active, list(df_Cas9$Day.c), FUN = mean)
mean.day.inactive_Cas9 <- aggregate(df_Cas9$p.inactive,list(df_Cas9$Day.c), FUN = mean)

mean.day.uncut_Cas9 <- aggregate(df_Cas9$p.uncut[df_Cas9$replicate == "Replicate2"], 
                                 list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate2"]), FUN = mean)
mean.day.active_Cas9 <- aggregate(df_Cas9$p.active[df_Cas9$replicate == "Replicate2"], 
                                  list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate2"]), FUN = mean)
mean.day.inactive_Cas9 <- aggregate(df_Cas9$p.inactive[df_Cas9$replicate == "Replicate2"], 
                                    list(df_Cas9$Day.c[df_Cas9$replicate == "Replicate2"]), FUN = mean)

df2_Cas9 <- data.frame(
  Day = rep(mean.day.uncut_Cas9[, 1], 3),
  value = c(mean.day.uncut_Cas9[, 2], mean.day.active_Cas9[, 2], mean.day.inactive_Cas9[, 2]),
  activity = c(rep("Original", nrow(mean.day.uncut_Cas9)), 
               rep("Active", nrow(mean.day.active_Cas9)), 
               rep("Inactive", nrow(mean.day.inactive_Cas9)))
)

df2_Cas9$activity <- factor(df2_Cas9$activity, levels = c("Original", "Active", "Inactive"))

ggplot(data = df2_Cas9, aes(x = as.numeric(as.vector(Day)), y = value, fill = activity)) +
  geom_area() +
  scale_fill_brewer(palette="Oranges")  +
  theme_classic() +
  labs(x = "Day", y = "Proportion of barcodes", fill = "Activity", title = "Cas9")+
  theme(plot.title = element_text(hjust = 0.5))


#     BE3:
# -------------
df_BE3 <- df[df$System == "BE3", ]
mean.day.uncut_BE3 <- aggregate(df_BE3$p.uncut[df_BE3$replicate == "Replicate2"], 
                                list(df_BE3$Day.c[df_BE3$replicate == "Replicate2"]), FUN = mean)
mean.day.active_BE3 <- aggregate(df_BE3$p.active[df_BE3$replicate == "Replicate2"], 
                                 list(df_BE3$Day.c[df_BE3$replicate == "Replicate2"]), FUN = mean)
mean.day.inactive_BE3 <- aggregate(df_BE3$p.inactive[df_BE3$replicate == "Replicate2"], 
                                   list(df_BE3$Day.c[df_BE3$replicate == "Replicate2"]), FUN = mean)
df2_BE3 <- data.frame(
  Day = rep(mean.day.uncut_BE3[, 1], 3),
  value = c(mean.day.uncut_BE3[, 2], mean.day.active_BE3[, 2], mean.day.inactive_BE3[, 2]),
  activity = c(rep("Original", nrow(mean.day.uncut_BE3)), 
               rep("Active", nrow(mean.day.active_BE3)), 
               rep("Inactive", nrow(mean.day.inactive_BE3)))
)
df2_BE3$activity <- factor(df2_BE3$activity, levels = c("Original", "Active", "Inactive"))
ggplot(data = df2_BE3, aes(x = as.numeric(as.vector(Day)), y = value, fill = activity)) +
  geom_area() +
  scale_fill_brewer(palette="BuGn") +
  theme_classic() +
  labs(x = "Day", y = "Proportion of barcodes", fill = "Activity", title = "BE3")+
  theme(plot.title = element_text(hjust = 0.5))


#   Figure 2d : . Barcode diversity over time
# =======================================================================================
ggplot(metadata, aes(x = factor(Day.c), y = Diversity, fill = System)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#2ca25f",  "#e6550d")) + 
  facet_wrap(~ System) +  
  labs(x = "Day", y = "Barcode Diversity", title = "") +
  theme_classic()


#   Supplementary Figure 4: Proportion of original barcodes along with barcode diversity over time
# =======================================================================================
data_r1 <- metadata %>%
  filter(replicate == "Replicate1", Coverage > 200)

# Combine Original Sequence and Barcode Diversity into a long format
data_combined <- data_r1 %>%
  pivot_longer(cols = c(p.uncut, Diversity), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(Metric = case_when(
    Metric == "p.uncut" ~ "% Original sequence",
    Metric == "Diversity" ~ "Barcode Diversity"
  ))

# Plot all metrics in the same plot
ggplot(data_combined, aes(x = as.factor(Day.c), y = Value, 
                          color = interaction(System, Metric), 
                          group = interaction(System, Metric))) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) +  
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9.% Original sequence" = "#e6550d", 
                                "BE3.% Original sequence" = "#2ca25f", 
                                "Cas9.Barcode Diversity" = "#f27059", 
                                "BE3.Barcode Diversity" = "#52796f")) +
  scale_shape_manual(values = c(16, 17)) +  
  scale_y_continuous(
    name = "% Original sequence", 
    sec.axis = sec_axis(~., name = "Barcode Diversity") # Segundo eje
  ) +
  labs(title = "",
       x = "Day") +
  facet_wrap(~gRNA) +  
  theme_minimal() +
  theme(legend.position = "right",  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.y.right = element_text(size = 14, color = "black"),  # Color del eje secundario
        strip.text = element_text(size = 12, face = "bold"))

#   Supplementary Figure 5: Average of indels per read grouped by system (BE3 and canonical Cas9)
# =======================================================================================
ggplot(data_r1, aes(x = System, y = Indel, fill = System)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  labs(title = "Proportion of indels in Cas9 vs BE3",
       x = "System",
       y = "% Indels") +
  theme_minimal()

wilcox.test(Indel ~ System, data = data_r1)
