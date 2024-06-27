setwd("D:/RStudio/Projects/SARCOMA/raw_data/Cibersort") # set wd

#get packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("DESeq2")
library("tidyverse")
library("readxl")
library("ggplot2")
library("pheatmap")
library("ComplexHeatmap")
library("EnhancedVolcano")
library("vegan")
library("org.Hs.eg.db")
library("dplyr")
library("heatmap.plus")
library("heatmap.plus")
library("circlize")
library("pheatmap")
library("dplyr")
library("reshape2")
library("ggsignif")

cell_proportions <- read.csv("CIBERSORTx_results.csv", row.names = 1, header=T) # get CIBERSORTx results
cell_proportions<- cell_proportions[,order(colnames(cell_proportions))]
cell_proportions <- cell_proportions %>%
  dplyr::select(-P.value, -Correlation, -RMSE)


col_MS<-read.csv("col_LMS.csv", header=T, sep = ";", row.names = 1) #clinical data, modified

col_LMS_subset <- dplyr::select(col_LMS, H2AX_status)
col_LMS_subset$Sample_ID <- rownames(col_LMS_subset)

# Add the rownames as a new column to the cell_proportions data frame for merging
cell_proportions$Sample_ID <- rownames(cell_proportions)

# Merge the data frames on the Sample_ID column
merged_data <- merge(col_LMS_subset, cell_proportions, by = "Sample_ID")


# Melt the data for ggplot2, including H2AX_status
melted_data <- melt(merged_data, id.vars = c("Sample_ID", "H2AX_status"))
# Rename the melted data columns for clarity
colnames(melted_data) <- c("Sample_ID", "H2AX_status", "Cell_Type", "Proportion")
melted_data$Cell_Type <- factor(melted_data$Cell_Type, levels = sort(unique(melted_data$Cell_Type)))

# Perform t-tests and collect p-values
cell_types <- colnames(cell_proportions)[-ncol(cell_proportions)]  # Exclude Sample_ID column

t_test_results <- data.frame(Cell_Type = character(), P_Value = numeric(), Significance = character())

for (cell_type in cell_types) {
  high_group <- merged_data %>% filter(H2AX_status == "1") %>% pull(cell_type)
  low_group <- merged_data %>% filter(H2AX_status == "0") %>% pull(cell_type)
  
  t_test <- t.test(high_group, low_group)
  
  significance <- ""
  if (t_test$p.value <= 0.001) {
    significance <- "***"
  } else if (t_test$p.value <= 0.01) {
    significance <- "**"
  } else if (t_test$p.value <= 0.05) {
    significance <- "*"
  }
  
  t_test_results <- rbind(t_test_results, data.frame(
    Cell_Type = cell_type,
    P_Value = t_test$p.value,
    Significance = significance
  ))
}

# Print t-test results
print(t_test_results)

# Convert H2AX_status to factor
melted_data$H2AX_status <- as.factor(melted_data$H2AX_status)

if (length(missing_levels) > 0) {
  cat("Warning: Some Cell_Type have missing levels of H2AX_status:\n")
  print(missing_levels)
}

# Plot the data
png(file="cibersort_barplot.png", width=22, height=12, units = "cm", res=1200)
ggplot(melted_data, aes(x = Cell_Type, y = Proportion, fill = H2AX_status)) +
  geom_boxplot(position = position_dodge(width = 1.0), orientation = "horizontal") +  # Adjust the width parameter here
  scale_fill_manual(values = c("0" = "blue", "1" = "darkorange"), name = "H2AX Status",
                    labels = c("Low", "High")) +
  labs(x = "Cell Type", y = "Relative Proportion", title = "Cell Type Proportions by H2AX Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

median(cell_proportions$Macrophages.M1)
# Ensure the rownames of both data frames are correctly set
rownames(cell_proportions) <- cell_proportions$rownames_column
rownames(col_LMS) <- col_LMS$rownames_column

# Select the "Macrophages.M1" column from cell_proportions
macrophages_m1 <- cell_proportions["Macrophages.M1"]

# Combine the data frames by rownames
col_LMS <- cbind(col_LMS, macrophages_m1[rownames(col_LMS),, drop = FALSE])

# If you need to remove any rownames column used in the process, uncomment the following lines:
# col_LMS$rownames_column <- NULL
# cell_proportions$rownames_column <- NULL


col_LMS$m1_status<-ifelse(col_LMS$Macrophages.M1 >=0.001, "1", "0")
col_LMS$survival_m <- col_LMS$survival_d/30
write.table(col_LMS, file = "col_LMS_Immunescores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

km_m1 <- survfit(Surv(survival_m, vital_status) ~m1_status,
                      data=col_LMS,
                      type = "kaplan-meier")

ggsurvplot(km_m1, data=col_LMS, conf.int = F, risk.table = F,
           pval = T, pval.method = F, palette = c("black", "purple3"), 
           legend = "right", 
           legend.title = "M1 macrophage proportions", 
           legend.lab = c("< median", "≥ median"),
           xlab = "Time [months]", ylab = "Survival probability",
           ggtheme = theme_survminer())
median(col_LMS$Macrophages.M1)

#let's draw a Heatmap
library(ComplexHeatmap)

# Transpose the cell proportions matrix
cp_heatmap<- dplyr::select(cell_proportions, -Sample_ID)
cp_heatmap<-t(cp_heatmap)
cp_heatmap<-data.frame(cp_heatmap)
col_LMS_heatmap<-col_LMS_subset

cp_heatmap_sort <- cp_heatmap[, intersect(rownames(col_LMS_heatmap), colnames(cp_heatmap))]
cp_heatmap<-cp_heatmap[,colnames(cp_heatmap_sort)]
col_LMS_heatmap<-col_LMS_heatmap[colnames(cp_heatmap),]

cp_heatmap<-cp_heatmap <- cp_heatmap[order(rownames(cp_heatmap)), ]
dim(cp_heatmap)

# Create a dataframe for column split
colsplit <- data.frame(H2AX_status = col_LMS_heatmap$H2AX_status)
rownames(colsplit) <- rownames(col_LMS_heatmap)



# Define colors for column split
levels_to_colors <- list(H2AX_status = c("0" = "blue", "1" = "darkorange"))

# Create HeatmapAnnotation object
ha <- HeatmapAnnotation(df = colsplit, col = levels_to_colors, na_col = "white")

# Create the heatmap
heat <- Heatmap(as.matrix(cp_heatmap), top_annotation = ha, column_split = colsplit,
                name = "Relative Proportion",
                cluster_columns = FALSE,  # Do not cluster columns (samples)
                cluster_rows = FALSE,      # Cluster rows (cell types)
                column_title = "Sample",
                row_title = "Cell Type",
                show_row_names = TRUE,    # Show row names (cell types)
                show_column_names = FALSE, 
                row_names_side = "left",
                # Do not show column names (sample IDs)
                col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.5), c("black","indianred", "indianred2", "red", "red4")))  # Set legend direction

print(heat)


library(dplyr)
library(reshape2)

# Load and process cell proportions data
cell_proportions <- read.csv("CIBERSORT_rel.csv", row.names = 1, header = TRUE)
cell_proportions <- cell_proportions[, order(colnames(cell_proportions))]
cell_proportions <- cell_proportions %>% dplyr::select(-P.value, -Correlation, -RMSE)

colnames(cell_proportions) <- c("B cells, memory", "B cells, naive", "Dendritic cells, activated", 
                                "Dendritic cells, resting", "Eosinophils", "M0 Macrophages", 
                                "M1 Macrophages", "M2 Macrophages", "Mast cells, activated", 
                                "Mast cells, resting", "Monocytes", "Neutrophils", 
                                "NK cells, activated", "NK cells, resting", "Plasma cells", 
                                "CD4 T memory cells, activated", "CD4 T memory cells, resting", 
                                "CD4 T cells, naive", "CD8 T cells", "T follicular helper cells", 
                                "T cells gamma delta", "T cells, regulatory")

# Load and process metadata
col_LMS <- read.csv("col_LMS.csv", header = TRUE, sep = ";", row.names = 1)
col_LMS_subset <- dplyr::select(col_LMS, H2AX_status)
col_LMS_subset$Sample_ID <- rownames(col_LMS_subset)

# Add Sample_ID column to cell_proportions for merging
cell_proportions$Sample_ID <- rownames(cell_proportions)

# Merge data frames on Sample_ID
merged_data <- merge(col_LMS_subset, cell_proportions, by = "Sample_ID")

# Melt the data for ggplot2, including H2AX_status
melted_data <- melt(merged_data, id.vars = c("Sample_ID", "H2AX_status"))
colnames(melted_data) <- c("Sample_ID", "H2AX_status", "Cell_Type", "Proportion")
melted_data$Cell_Type <- factor(melted_data$Cell_Type, levels = sort(unique(melted_data$Cell_Type)))

# Perform t-tests and collect p-values along with means for each subgroup
cell_types <- colnames(cell_proportions)[-ncol(cell_proportions)]  # Exclude Sample_ID column

t_test_results <- data.frame(Cell_Type = character(), 
                             P_Value = numeric(), 
                             Significance = character(), 
                             Mean_H2AX_low = numeric(), 
                             Mean_H2AX_high = numeric())

for (cell_type in cell_types) {
  high_group <- merged_data %>% filter(H2AX_status == "1") %>% pull(cell_type)
  low_group <- merged_data %>% filter(H2AX_status == "0") %>% pull(cell_type)
  
  t_test <- t.test(high_group, low_group)
  
  significance <- ""
  if (t_test$p.value <= 0.001) {
    significance <- "***"
  } else if (t_test$p.value <= 0.01) {
    significance <- "**"
  } else if (t_test$p.value <= 0.05) {
    significance <- "*"
  }
  
  t_test_results <- rbind(t_test_results, data.frame(
    Cell_Type = cell_type,
    P_Value = t_test$p.value,
    Significance = significance,
    Mean_H2AX_low = mean(low_group, na.rm = TRUE),
    Mean_H2AX_high = mean(high_group, na.rm = TRUE)
  ))
}

# Print t-test results
print(t_test_results)

# Specifically check for Eosinophils and resting Dendritic cells
eosinophils_result <- t_test_results %>% filter(Cell_Type == "Eosinophils")
dendritic_resting_result <- t_test_results %>% filter(Cell_Type == "Dendritic cells, resting")

print(eosinophils_result)
print(dendritic_resting_result)
