# Load the PBMC dataset
pbmc4.data <- Read10X(data.dir = "data/BCP004")


# Initialize the Seurat object with the raw (non-normalized data).
pbmc4 <- CreateSeuratObject(counts = pbmc4.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc4[["percent.mt"]] <- PercentageFeatureSet(pbmc4, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#BCP004
pbmc4 <- subset(pbmc4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)


pbmc4 <- NormalizeData(pbmc4)

pbmc4 <- FindVariableFeatures(pbmc4, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc4), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes4 <- rownames(pbmc4)
pbmc4 <- ScaleData(pbmc4, features = all.genes4)

pbmc4 <- RunPCA(pbmc4, features = VariableFeatures(object = pbmc4))

# Examine and visualize PCA results a few different ways
print(pbmc4[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc4, dims = 1:2, reduction = "pca")

DimPlot(pbmc4, reduction = "pca") + NoLegend()

DimHeatmap(pbmc4, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc4, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc4)

pbmc4 <- FindNeighbors(pbmc4, dims = 1:10)

#change the resolution of the clusters to group more/less cells together
pbmc4 <- FindClusters(pbmc4, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc4), 5)

pbmc4 <- RunUMAP(pbmc4, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc4, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc4, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc4, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc4.markers <- FindAllMarkers(pbmc4, only.pos = TRUE)
pbmc4.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc4, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)



