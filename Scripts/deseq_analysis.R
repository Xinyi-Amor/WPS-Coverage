library(DESeq2)
library(ggplot2)

coverage <- read.delim("output_processed.txt", header=FALSE, row.names=1, sep="\t")
# coverage <- coverage[which(rowSums(coverage)>1)]
condition <- factor(c("C","C","C","S","S","S"))

coldata <- data.frame(row.names=colnames(coverage),condition)

dds <- DESeqDataSetFromMatrix(countData = coverage, colData = coldata, design = ~condition)
dds <- DESeq(dds)
#vsdata <- vst(dds, blind=FALSE)

#pca_plot <- plotPCA(vsdata, intgroup = "condition")
#png("PCA_plot.png")  # Change the file format if needed
#print(pca_plot)
#dev.off()

res <- results(dds, contrast = c("condition", "S", "C"))
results <- data.frame(
  gene_name = rownames(res),
  logFC = res$log2FoldChange,
  p_value = res$pvalue
)

write.table(results, file = "deseq_coverage_processed_logFC.txt", sep = "\t", row.names = FALSE)

# Plot log-ratio vs. average expression
#plotdes <- plot(res$log2FoldChange, -log10(res$pvalue), xlim = c(-2, 2), xlab = "Log2 Fold Change", ylab = "-log10(P-value)")
# Add labels with row names
#text(res$log2FoldChange, -log10(res$pvalue), labels = rownames(res), pos = 4, offset = 0.5, cex = 0.2)

#png("des_plot.png")
#print(plotdes)
#dev.off()

#library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
#library(RColorBrewer) # for a colourful plot
#library("org.Hs.eg.db", character.only = TRUE)
#library(msigdbr)
#library(fgsea)
#m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

#fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#rankings <- sign(res$log2FoldChange) * (-log10(res$pvalue))
#names(rankings) <- rownames(res)
#rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
#GSEAres <- fgsea(fgsea_sets, stats = rankings) #ranks1 comes from your differential gene expression data
#head(GSEAres[order(pval), ])
#plot <- plotEnrichment(fgsea_sets[["REACTOME_PTEN_REGULATION"]],
#               rankings) + labs(title="PTEN Regulation")
#ggsave("output.png", plot, width = 8, height = 6)
#topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
#topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
#topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
#plot_gsea <- plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#png("gsea_plot.png")
#print(plot_gsea)
#dev.off()
