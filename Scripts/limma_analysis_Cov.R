library(limma)
library(ggplot2)
library(edgeR)

# Step 1: Read in the data
count_data <- read.delim("output_coverage.txt", header = FALSE, row.names = 1, sep = "\t")
# Step 2: Normalize the data     
d <- DGEList(counts = count_data)
v <- voom(d, normalize.method = "scale")

count_data <- v$E

# Step 2: Define groups
group1_cols <- 1:3  # Columns corresponding to group 1
group2_cols <- 4:6  # Columns corresponding to group 2

# Step 3: Fit a linear model
design <- model.matrix(~0 + factor(c(rep("Group1", length(group1_cols)), rep("Group2", length(group2_cols)))), data = data.frame(count_data))
colnames(design) <- c("Group1", "Group2")
fit <- lmFit(count_data, design)

# Step 4: Test for differential expression
contrast.matrix <- makeContrasts(Group2 - Group1, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

# Step 5: Extract and process results
logFC <- fit$coefficients[, "Group2 - Group1"]
p_values <- fit$p.value[, "Group2 - Group1"]
results <- data.frame(
  gene_name = rownames(count_data),
  logFC = logFC,
  p_value = p_values
)

write.table(results, file = "limma_coverage_logFC_P_Scal.txt", sep = "\t")

#library("org.Hs.eg.db", character.only = TRUE)
#library(msigdbr)
#library(fgsea)

#m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

#fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#rankings <- sign(results$logFC)*(-log10(results$p_value))

# Extract part before the period in gene names
#gene_names <- sub("\\..*", "", results$gene_name)

#names(rankings) <- gene_names

#rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
#GSEAres <- fgsea(fgsea_sets, stats = rankings) #rankings comes from your differential gene expression data
#head(GSEAres[order(pval), ])

# Create a volcano plot
#volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(p_value))) + geom_point(alpha = 0.5) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = "Coverage Plot")
#volcano_plot <- volcano_plot + geom_text(aes(label = rownames(results)), , vjust = -0.5, size = 1.5)
#ggsave("new_coverage_plot_Scal.png", plot = volcano_plot, width = 6, height = 4, units = "in", dpi = 300)

