#Gene signatures in PCA clusters Color Map

library(Seurat)
library(ggplot2)

# Load the dataset
idhwt_data <- read.csv("C:/Users/uabic/Desktop/IDHWT.csv", row.names = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = idhwt_data)

# Log normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Plot PCA clustering
pca_plot <- DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Clustering")

# Define the genes of interest
genes_of_interest <- c("AQP4", "S100B", "GFAP", "SOX2", "SOX11", "PROM1", "PDGFRA", "EGFR", "CD44", "CD74", "NES", "TMEM119")

library (patchwork)
# Generate expression plots
expression_plots <- lapply(genes_of_interest, function(gene) {
  FeaturePlot(seurat_obj, features = gene, reduction = "pca") + ggtitle(paste("Expression of", gene))
})

# Combine plots into two panels
combined_plot_1 <- wrap_plots(expression_plots[1:6])
combined_plot_2 <- wrap_plots(expression_plots[7:12])


# Save the PCA plot
jpeg("PCA_Clustering.jpeg", width = 800, height = 800, quality = 100)
print(pca_plot)
dev.off()

# Save the combined expression plots
jpeg("Gene_Expression_PCA_Panel1.jpeg", width = 1920, height = 1080, quality = 100)
print(combined_plot_1)
dev.off()

jpeg("Gene_Expression_PCA_Panel2.jpeg", width = 1920, height = 1080, quality = 100)
print(combined_plot_2)
dev.off()
