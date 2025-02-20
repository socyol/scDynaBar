##############################################
###        Bulk analysis Figures          ###
##############################################

#_______________________________________
# 1. Load libraries and scripts
# ======================================
source("scripts/settings.R")

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

#   Figure 1c :  . Proportion of barcodes based on their mutational profile: active, inactive); and original 
# =======================================================================================

# df$p.active.uncut.PAM <- df$active.UncutPAM.qualiy/df$Coverage
# df$p.inactive.cut  <- df$inactive.PAM.quality/df$Coverage

df <-  df[df$Coverage > 200,]

aggregate(df$p.uncut.quality[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.uncut
aggregate(df$p.active.uncut.PAM[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.active.uncutPAM
aggregate(df$p.inactive.cut[df$replicate == "Replicate1"], list(df$Day.c[df$replicate == "Replicate1"]), FUN=mean) -> mean.day.inactive
df2 <- data.frame( Day= rep(mean.day.inactive[,1], 3), value= c( mean.day.inactive[,2],mean.day.active.uncutPAM[,2],mean.day.uncut[,2] ),
                   activity= c( rep("inactive", 6), rep( "activeUncutPAM", 6), rep( "uncut", 6)))
act <- factor(df2$activity, levels=c(  "uncut", "activeUncutPAM", "inactive" ))
ggplot(data=df2, aes(x=as.numeric(as.vector(Day)), y=value, fill= act)) +
  geom_area() +   scale_fill_brewer(palette="Blues") + theme_classic()


#   Figure 1d :  Proportion of original barcodes over time across different gRNAs
# =======================================================================================
ggplot(data_r1, aes(x = as.factor(Day.c), y = p.uncut, color = System, group = System)) +
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
ggplot(data_r1, aes(x = as.factor(Day.c), y = Diversity.y, color = System, group = System)) +
  geom_point(aes(shape = System), size = 3, position = position_dodge(width = 0.5)) + 
  geom_line(size = 1.2) + 
  scale_color_manual(values = c("Cas9" = "#e6550d", "BE3" = "#2ca25f")) + 
  scale_shape_manual(values = c(16, 17)) +  
  labs(title = "",
       x = "",
       y = "Barcode Diversity") +
  facet_wrap(~gRNA) +  
  theme_minimal() +
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
    gRNA %in% c("g6", "g7") ~ "26bp",
    TRUE ~ NA_character_  
  ))
ggplot(data_r1_bp, aes(x = factor(Day.c), y = p.uncut)) +
  geom_boxplot() +
  facet_wrap(~ spacer) +  
  scale_fill_brewer(palette = "Purples",  direction = 1) +
  theme_classic()+
  labs(x = "Day", y = "% Original sequence", title = "")

#   Figure 1g : Percentage of the original nucleotide over time, aligning the spacer relative to the PAM sequence (Cas9 system)
# =======================================================================================
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

tmp2 <- tmp[df$Experiment =="Cas9" ,]
tmp3 <- rbind( colMeans(tmp2[ rownames(tmp2) == 0, ]), colMeans(tmp2[ rownames(tmp2) == 2, ] ), colMeans(tmp2[ rownames(tmp2) == 4, ]),
               colMeans(tmp2[ rownames(tmp2) == 8, ]), colMeans(tmp2[ rownames(tmp2) == 10, ]), colMeans(tmp2[ rownames(tmp2) == 18, ]),
               colMeans(tmp2[ rownames(tmp2) == 25, ]), colMeans(tmp2[ rownames(tmp2) == 32, ]))


rampCol2 <- colorRampPalette(c( "red", "#f03b20", "#feb24c" , "#ffffcc","white"))(n = 1000)
library(gplots)

heatmap.2( tmp3 , na.color ="#f1a34", Colv = FALSE, Rowv =FALSE ,trace="none" , col=viridis )


#   Figure 1h : Percentage of C>T mutations over time
# =======================================================================================
ggplot(df[ df$replicate == "Replicate1",], aes(x=Day.c, y=C.T,  fill=Experiment)) + geom_boxplot( alpha=0.9) +
  theme_classic()  +    ylim(0, 0.1) +
  scale_fill_manual(values=c("#2ca25f",  "#e6550d"))






#   Figure 2b :  Pearson’s correlation coefficient for reproducibility between replicates from days 0 to 18.
# =======================================================================================
## correlation plot
replicate1 <- c("8227", "8238")
replicate2<- c( "8467", "8468" )

df1 <- df[df$replicate == "Replicate1" &  df$Reads_QC_amplicon > 200,]
df2 <- df[df$replicate == "Replicate2" & df$Reads_QC_amplicon > 200,]

df1 <- df[df$replicate == "Replicate1" &  df$Coverage > 200,]
df2 <- df[df$replicate == "Replicate2" & df$Coverage > 200,]

data.cor <-c()
data.mean  <-c()
for ( i in  c(1:nrow(df1))) {
  day= df1$Day.c[i]
  gRNA= df1$gRNA[i]
  exp= df1$Experiment[i]
  value <- which(df2$Day.c == day & df2$gRNA== gRNA & df2$Experiment == exp)  
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


# Número total de combinaciones únicas de gRNA, System y Day
num_combinations <- n_distinct(data_table[, c("Day", "gRNA", "System")]) # 33

# Número total de días y replicates considerados
num_days <- length(unique(data_table$Day)) # 6
num_replicates <- 2  # Ya que son replicate 1 y replicate 2

cat(" Total unique combinations (gRNA-System-Day):", num_combinations, "\n",
    "Total days analyzed:", num_days, "\n",
    "Total replicates per combination:", num_replicates, "\n")

write.csv(data_table, "TFM/paper_REVIEWER-COMENTS/replicates_bulk_table.csv", row.names = FALSE)


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

#     Cas9:
# -------------
aggregate(df$p.uncut.quality[df$Experiment == "Cas9" & df$replicate == "Replicate1"], list(df$Day.c[df$Experiment == "Cas9"& df$replicate == "Replicate1" ]), FUN=mean) -> mean.day.uncut
aggregate(df$p.active.uncut.PAM[df$Experiment == "Cas9"  & df$replicate == "Replicate1"], list(df$Day.c[df$Experiment == "Cas9"  & df$replicate == "Replicate1"]), FUN=mean) -> mean.day.active.uncutPAM
aggregate(df$p.inactive.cut[df$Experiment == "Cas9" & df$replicate == "Replicate1"], list(df$Day.c[df$Experiment == "Cas9" & df$replicate == "Replicate1"]), FUN=mean) -> mean.day.inactive
df2 <- data.frame( Day= rep(mean.day.inactive[,1], 3), value= c( mean.day.inactive[,2],mean.day.active.uncutPAM[,2],mean.day.uncut[,2] ),
                   activity= c( rep("inactive", 6), rep( "activeUncutPAM", 6), rep( "uncut", 6)))

act <- factor(df2$activity, levels=c("uncut" , "activeUncutPAM",  "inactive" ))
ggplot(data=df2, aes(x=as.numeric(as.vector(Day)), y=value, fill= act)) +
  geom_area() +   scale_fill_brewer(palette="Oranges") + theme_classic()

#     BE3:
# -------------
aggregate(df$p.uncut.quality[df$replicate == "Replicate2" & df$Experiment == "BE3"], list(df$Day.c[df$Experiment == "BE3" & df$replicate == "Replicate2" ]), FUN=mean) -> mean.day.uncut
aggregate(df$p.active.uncut.PAM[df$replicate == "Replicate2" & df$Experiment == "BE3"], list(df$Day.c[df$Experiment == "BE3" & df$replicate == "Replicate2" ]), FUN=mean) -> mean.day.active.uncutPAM
aggregate(df$p.inactive.cut[df$replicate == "Replicate2" & df$Experiment == "BE3"], list(df$Day.c[df$Experiment == "BE3" & df$replicate == "Replicate2" ]), FUN=mean) -> mean.day.inactive
df2 <- data.frame( Day= rep(mean.day.inactive[,1], 3), value= c( mean.day.inactive[,2],mean.day.active.uncutPAM[,2],mean.day.uncut[,2] ),
                   activity= c( rep("inactive", 7), rep( "activeUncutPAM", 7), rep( "uncut", 7)))
act <- factor(df2$activity, levels=c(    "uncut", "activeUncutPAM", "inactive"))
ggplot(data=df2, aes(x=as.numeric(as.vector(Day)), y=value, fill= act)) +
  geom_area() +   scale_fill_brewer(palette="BuGn") + theme_classic()


#   Figure 2d : . Barcode diversity over time
# =======================================================================================
ggplot(data_r1, aes(x = factor(Day.c), y = Diversity.y, fill = System)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#2ca25f",  "#e6550d")) + 
  facet_wrap(~ System) +  
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

