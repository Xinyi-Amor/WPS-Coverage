library(ggplot2)
# Step 1: Read data from text files
x_values <- read.table("deseq_coverage_processed_logFC.txt", header = TRUE, sep ='\t') 
y_values <- read.table("limma_coverage_processed_logFC.txt", header = TRUE, sep ='\t')

# Step 2: Extract columns
x <- x_values$logFC  # Assuming the column you want to use from x_values is the second column
y <- y_values$logFC  # Assuming the column you want to use from y_values is the second column

# Step 3: Create a data frame
data <- data.frame(
        gene_name = x_values$gene_name,
        deseq_logFC = x,
        limma_logFC = y
)  

# Step 4: Create the plot using ggplot
plot <- ggplot(data, aes(x = deseq_logFC, y = limma_logFC)) +
  geom_point(color = "blue", size = 3) +
  labs(x = "deseq coverage logFC", y = "limma coverage logFC", title = "Scatter Plot")
plot <- plot + geom_text(aes(label = gene_name), vjust = -0.5, size = 1.5)
# Step 5: Save the plot as an image file
ggsave("limma_deseq_plot.png", plot, width = 6, height = 6, dpi = 300)
