count_data <- read.delim("processed.txt", header = FALSE, row.names = 1, sep="\t")
print(head(count_data))
group1_cols <- 1:3  # Columns corresponding to group 1
group2_cols <- 4:6  # Columns corresponding to group 2

# Calculate mean expression for each group
mean_expression_group1 <- rowMeans(count_data[, group1_cols])
mean_expression_group2 <- rowMeans(count_data[, group2_cols])
# Compute log fold change
logFC <- log2(mean_expression_group2 / mean_expression_group1)
print(head(logFC))
# Initialize vectors to store p-values
gene_names <- rownames(count_data)
p_values <- numeric(length(gene_names))

# Perform t-test for each gene
for (i in seq_along(gene_names)) {
  group1_values <- count_data[i, group1_cols]
  group2_values <- count_data[i, group2_cols]
  t_test_result <- t.test(group2_values, group1_values)
  p_values[i] <- t_test_result$p.value
}

# Create a data frame with gene names, logFC, and p-values
results <- data.frame(
  gene_name = gene_names,
  logFC = logFC,
  p_value = p_values
)

# Write the data frame to a tab-delimited text file (e.g., output.txt)
write.table(results, file = "test2_coverage_logFC_P.txt", sep = "\t")

library(ggplot2)
# Create a volcano plot using the results data frame
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(p_value))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "Log2 Fold Change", y = "-log10(P_value)", title = "Coverage Plot")
volcano_plot <- volcano_plot + geom_text(aes(label = results$gene_name), vjust = -0.5, size = 1.5)
ggsave("test2.png", plot = volcano_plot, width = 6, height = 4, units = "in", dpi = 300)


