library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc3.data <- Read10X(data.dir = "/data/BCP003")


# Initialize the Seurat object with the raw (non-normalized data).
pbmc3 <- CreateSeuratObject(counts = pbmc3.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc3[["percent.mt"]] <- PercentageFeatureSet(pbmc3, pattern = "^MT-")
  
# Visualize QC metrics as a violin plot
VlnPlot(pbmc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
plot1 <- FeatureScatter(pbmc3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#BCP003
pbmc3 <- subset(pbmc3, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 5)


pbmc3 <- NormalizeData(pbmc3)

pbmc3 <- FindVariableFeatures(pbmc3, selection.method = "vst", nfeatures = 6500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc3), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes3 <- rownames(pbmc3)
pbmc3 <- ScaleData(pbmc3, features = all.genes3)

pbmc3 <- RunPCA(pbmc3, features = VariableFeatures(object = pbmc3))

# Examine and visualize PCA results a few different ways
print(pbmc3[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc3, dims = 1:2, reduction = "pca")

DimPlot(pbmc3, reduction = "pca") + NoLegend()

DimHeatmap(pbmc3, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc3, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc3)

pbmc3 <- FindNeighbors(pbmc3, dims = 1:10)

#change the resolution of the clusters to group more/less cells together
pbmc3 <- FindClusters(pbmc3, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc3), 5)

pbmc3 <- RunUMAP(pbmc3, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc3, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc3, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc3, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc3.markers <- FindAllMarkers(pbmc3, only.pos = TRUE)
pbmc3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc3, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)



