# Data Testing and Cleaning 
data1 <- read.table("Downloads/Brain_Gene_Seq/csv/amygdala_reads.csv", header = TRUE, sep = "\t")
data2 <- read.table("Downloads/Brain_Gene_Seq/csv/hippocampus_reads.csv", header = TRUE, sep = "\t")
print(head(data1))

#Sets Threshold for datset - Expressed Genes must be > 150000 Total Count across Samples
filter_low_expressed_genes <- function(input_file, output_file, threshold = 150000, skip_lines = 0, row_names_col = 1) {
  data <- read.csv(input_file, skip = skip_lines, row.names = row_names_col, check.names = FALSE)
  numeric_cols <- sapply(data, is.numeric)
  if (sum(numeric_cols) == 0) {
    metadata_cols <- 1
    numeric_data <- as.data.frame(lapply(data[ , -(metadata_cols), drop=FALSE], as.numeric), row.names = rownames(data))
  } else {
    metadata_cols <- which(!numeric_cols)
    numeric_data <- data[ , numeric_cols, drop=FALSE]
  }
  gene_sums <- rowSums(numeric_data, na.rm = TRUE)
  keep_rows <- gene_sums >= threshold
  filtered_numeric <- numeric_data[keep_rows, , drop=FALSE]
  filtered_metadata <- if (length(metadata_cols) > 0) data[keep_rows, metadata_cols, drop=FALSE] else NULL
  filtered_data <- cbind(filtered_metadata, filtered_numeric)
  write.csv(filtered_data, output_file, quote = FALSE, row.names = TRUE)
  cat(sum(!keep_rows), "rows removed; remaining rows:", sum(keep_rows), "\n")
  cat("Filtered dataset saved to", output_file, "\n")
}

#Reformats csv files to guarantee elements are comma delimited
raw_lines <- readLines("top20_amygdala_genes.csv")
split_lines <- strsplit(raw_lines, ",")
data <- do.call(rbind, lapply(split_lines, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
write.csv(data, "top20_amygdala_genes.csv", row.names = FALSE, quote = FALSE)
data1 <- read.table("top20_amygdala_genes.csv", header = TRUE)
print(head(data1))

# CLEANING

filter_low_expressed_genes("Downloads/Brain_Gene_Seq/csv/amygdala_reads.csv", "amygdala_clean.csv")
filter_low_expressed_genes("Downloads/Brain_Gene_Seq/csv/hippocampus_reads.csv", "hippocampus_clean.csv")
data1 <- read.table("Downloads/Brain_Gene_Seq/csv/amygdala_reads.csv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data2 <- read.table("Downloads/Brain_Gene_Seq/csv/hippocampus_reads.csv", header = TRUE, sep = "\t")

data <- read.csv("hippocampus_clean_formatted.csv", row.names = 1)
data_num <- as.data.frame(lapply(data, as.numeric), row.names = rownames(data))
gene_sums <- rowSums(data_num, na.rm = TRUE)
print(head(gene_sums))

amygdala <- read.csv("amygdala_clean_formatted.csv", row.names = 1, check.names = FALSE)
hippocampus <- read.csv("hippocampus_clean_formatted.csv", row.names = 1, check.names = FALSE)

shared_genes <- intersect(rownames(amygdala), rownames(hippocampus))
amygdala_shared <- amygdala[shared_genes, ]
hippocampus_shared <- hippocampus[shared_genes, ]
combined <- cbind(amygdala_shared, hippocampus_shared)
write.csv(combined, "combined_shared_genes.csv")

amygdala_data <- read.csv("amygdala_clean_formatted.csv", row.names = 1, check.names = FALSE)
amygdala_num <- as.data.frame(lapply(amygdala_data, as.numeric), row.names = rownames(amygdala_data))
amygdala_sums <- rowSums(amygdala_num, na.rm = TRUE)

# Get top 20 expressed genes in Amygdala
top20_amygdala <- sort(amygdala_sums, decreasing = TRUE)[1:20]
cat("Top 20 most expressed genes in Amygdala:\n")
print(top20_amygdala)

hippocampus_data <- read.csv("hippocampus_clean_formatted.csv", row.names = 1, check.names = FALSE)
hippocampus_num <- as.data.frame(lapply(hippocampus_data, as.numeric), row.names = rownames(hippocampus_data))
hippocampus_sums <- rowSums(hippocampus_num, na.rm = TRUE)

# Get top 20 expressed genes in Hippocampus
top20_hippocampus <- sort(hippocampus_sums, decreasing = TRUE)[1:20]
cat("\nTop 20 most expressed genes in Hippocampus:\n")
print(top20_hippocampus)

write.csv(data.frame(GeneID = names(top20_amygdala), TotalCount = top20_amygdala),
          "top20_amygdala_genes.csv", row.names = FALSE)

write.csv(data.frame(GeneID = names(top20_hippocampus), TotalCount = top20_hippocampus),
          "top20_hippocampus_genes.csv", row.names = FALSE)

library(ggplot2)
library(reshape2)

data <- read.csv("top20_amygdala_genes.csv", check.names = FALSE)

data_long <- melt(data, id.vars = "Tissue", variable.name = "Gene", value.name = "Expression")
genes <- unique(data_long$Gene)
for (gene_name in genes) {
  
  gene_data <- subset(data_long, Gene == gene_name)
  p <- ggplot(gene_data, aes(x = Tissue, y = Expression, fill = Tissue)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Expression of", gene_name), y = "Expression Count") +
    theme_minimal()
    print(p)
  
  #ggsave(filename = paste0(gene_name, "_barplot.png"), plot = p, width = 6, height = 4)
}

ggplot(data_long, aes(x = Gene, y = Expression, color = Tissue)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Dotplot of Top 20 Genes by Tissue")

library(pheatmap)

data_mat <- as.matrix(data[,-1])
rownames(data_mat) <- data$Tissue

pheatmap(data_mat, cluster_rows = FALSE, cluster_cols = TRUE,
         scale = "none", main = "Expression Heatmap: Amygdala vs Hippocampus")

library(ggplot2)

amygdala_data <- read.csv("top20_amygdala_genes.csv")

# Line plot for Amygdala
ggplot(amygdala_data, aes(x = reorder(GeneID, -TotalCount), y = TotalCount, group = 1)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 20 Expressed Genes in Amygdala", x = "Gene ID", y = "Total Expression Count")

hippocampus_data <- read.csv("top20_hippocampus_genes.csv")

# Line plot for Hippocampus
ggplot(hippocampus_data, aes(x = reorder(GeneID, -TotalCount), y = TotalCount, group = 1)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 20 Expressed Genes in Hippocampus", x = "Gene ID", y = "Total Expression Count")

top5_amygdala <- sort(amygdala_gene_means, decreasing = TRUE)[1:5]
print(top5_amygdala)
top5_hippocampus <- sort(hippocampus_gene_means, decreasing = TRUE)[1:5]
print(top5_hippocampus)


amygdala_data <- read.csv("amygdala_clean_formatted.csv", row.names = 1, check.names = FALSE)
amygdala_numeric <- as.data.frame(lapply(amygdala_data, as.numeric), row.names = rownames(amygdala_data))
amygdala_gene_means <- rowMeans(amygdala_numeric, na.rm = TRUE)
print(head(amygdala_gene_means))

hippocampus_data <- read.csv("hippocampus_clean_formatted.csv", row.names = 1, check.names = FALSE)
hippocampus_numeric <- as.data.frame(lapply(hippocampus_data, as.numeric), row.names = rownames(hippocampus_data))
hippocampus_gene_means <- rowMeans(hippocampus_numeric, na.rm = TRUE)
print(head(hippocampus_gene_means))

amy <- read.csv("amygdala_clean_formatted.csv", row.names = 1, check.names = FALSE)
hip <- read.csv("hippocampus_clean_formatted.csv", row.names = 1, check.names = FALSE)
amy_t <- as.data.frame(t(amy))
hip_t <- as.data.frame(t(hip))
common_genes <- intersect(colnames(amy_t), colnames(hip_t))
amy_t_common <- amy_t[, common_genes]
hip_t_common <- hip_t[, common_genes]
amy_t_common$Tissue <- "Amygdala"
hip_t_common$Tissue <- "Hippocampus"
df <- rbind(amy_t_common, hip_t_common)
df$Tissue <- factor(df$Tissue)

predictor_genes <- common_genes[1:3]
formula <- as.formula(paste("Tissue ~", paste(predictor_genes, collapse = " + ")))
model <- glm(formula, data = df, family = "binomial")
summary(model)
predicted_probs <- predict(model, type = "response")
predicted_labels <- ifelse(predicted_probs > 0.5, "Hippocampus", "Amygdala")
table(Predicted = predicted_labels, Actual = df$Tissue)

df_numeric <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
rownames(df_numeric) <- rownames(df)
df<- df_numeric[complete.cases(df_numeric), ]

tissue_labels <- df$Tissue
df$Tissue <- NULL 
dist_matrix <- dist(df)
hc <- hclust(dist_matrix, method = "complete")
plot(hc, labels = tissue_labels, main = "Sample clustering dendrogram")

