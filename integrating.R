# Load datasets
library(Seurat)

# Identify anchors
anchors <- FindIntegrationAnchors(object.list = list(pbmc3, pbmc4, pbmc5))

# Integrate datasets
integrated <- IntegrateData(anchorset = anchors)

# run standard anlaysis workflow
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated)

all.genesIntegrated <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genesIntegrated)

integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated), assay = "RNA")


integrated <- FindNeighbors(integrated, dims = 1:30, reduction = "pca")
integrated <- FindClusters(integrated, resolution = 0.1)

integrated <- IntegrateLayers(object = integrated, method = CCAIntegration, assay = "RNA", orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.1)

integrated <- RunUMAP(integrated, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(integrated, reduction = "umap")
