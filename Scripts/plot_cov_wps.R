library(ggplot2)
# Step 1: Read data from text files
x_values <- read.table("limma_wps_logFC_P.txt", header = TRUE, sep ='\t') 
y_values <- read.table("limma_coverage_logFC_P_Scal.txt", header = TRUE, sep ='\t')

# Step 2: Extract columns
x <- x_values$logFC  # Assuming the column you want to use from x_values is the second column
y <- y_values$logFC  # Assuming the column you want to use from y_values is the second column
#x <- x_values$p_value
#y <- y_values$p_value

# Step 3: Create a data frame
data <- data.frame(
        gene_name = x_values$gene_name,
        wps_p = x,
        cov_p = y
)  

# Step 4: Create the plot using ggplot
plot <- ggplot(data, aes(x = wps_p, y = cov_p)) +
  geom_point(color = "blue", size = 3) +
  labs(x = "limma wps logFC", y = "limma coverage logFC", title = "Scatter Plot")
plot <- plot + geom_text(aes(label = gene_name), vjust = -0.5, size = 1.5)
# Step 5: Save the plot as an image file
ggsave("limma_cov_wpsNoScal_plot_logFC.png", plot, width = 6, height = 6, dpi = 300)
