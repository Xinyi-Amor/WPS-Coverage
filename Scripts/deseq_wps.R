library(DESeq2)
library(ggplot2)

wps <- read.delim("wps_output.txt", header=FALSE, row.names=1, sep="\t")
# coverage <- coverage[which(rowSums(coverage)>1)]
condition <- factor(c("C","C","C","S","S","S"))

coldata <- data.frame(row.names=colnames(wps),condition)

dds <- DESeqDataSetFromMatrix(countData = wps, colData = coldata, design = ~condition)
dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)

pca_plot <- plotPCA(vsdata, intgroup = "condition")
png("PCA_plot.png")  # Change the file format if needed
print(pca_plot)
dev.off()

res <- results(dds, contrast = c("condition", "S", "C"))
# Plot log-ratio vs. average expression
plotdes <- plot(res$log2FoldChange, -log10(res$pvalue), xlim = c(-2, 2), xlab = "Log2 Fold Change", ylab = "-log10(P-value)")
# Add labels with row names
text(res$log2FoldChange, -log10(res$pvalue), labels = rownames(res), pos = 4, offset = 0.5, cex = 0.2)

png("des_plot.png")
print(plotdes)
dev.off()


