# Load the PBMC dataset
pbmc5.data <- Read10X(data.dir = "/data/BCP005")


# Initialize the Seurat object with the raw (non-normalized data).
pbmc5 <- CreateSeuratObject(counts = pbmc5.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc5[["percent.mt"]] <- PercentageFeatureSet(pbmc5, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#BCP005
pbmc5 <- subset(pbmc5, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)


pbmc5 <- NormalizeData(pbmc5)
pbmc5 <- NormalizeData(pbmc5, normalization.method = "LogNormalize", scale.factor = 100000)

pbmc5 <- FindVariableFeatures(pbmc5, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc5), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes5 <- rownames(pbmc5)
pbmc5 <- ScaleData(pbmc5, features = all.genes5)

pbmc5 <- RunPCA(pbmc5, features = VariableFeatures(object = pbmc5))

# Examine and visualize PCA results a few different ways
print(pbmc5[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc5, dims = 1:2, reduction = "pca")

DimPlot(pbmc5, reduction = "pca") + NoLegend()

DimHeatmap(pbmc5, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc5, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc5)

pbmc5 <- FindNeighbors(pbmc5, dims = 1:10)

#change the resolution of the clusters to group more/less cells together
pbmc5 <- FindClusters(pbmc5, resolution = 0.1)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc5), 5)

pbmc5 <- RunUMAP(pbmc5, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc5, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc5, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc5, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc5.markers <- FindAllMarkers(pbmc5, only.pos = TRUE)
pbmc5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc5, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)



