# Example pipeline for preprocessing and analyzing single-cell RNA-seq data using Seurat
#
# Assumes a SingleCellExperiment object 'sce' already processed through basic QC
# Steps include: creating a Seurat object, additional QC filtering, normalization,
# feature selection, dimensionality reduction, clustering, and cell type annotation.

library(Seurat)
library(SingleCellExperiment)

# Convert SingleCellExperiment to Seurat
seurat <- CreateSeuratObject(counts = counts(sce),
                             meta.data = as.data.frame(colData(sce)))

# Example QC: calculate mitochondrial gene percentage and filter cells
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
# Visualize QC metrics (optional)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Filter cells based on QC thresholds
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

# Identify highly variable genes
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

# Scale the data and run PCA
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- RunPCA(seurat, features = VariableFeatures(seurat))
ElbowPlot(seurat, ndims = 50)

# Choose number of PCs for downstream analyses
PCs <- 10

# Clustering and dimensionality reduction
seurat <- FindNeighbors(seurat, dims = 1:PCs)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:PCs)
seurat <- RunTSNE(seurat, dims = 1:PCs)

DimPlot(seurat, label = TRUE)

# Example canonical marker genes for annotation
marker.genes <- list(
  T.cell = c("CD3D", "CD3E", "TRAC"),
  Monocyte = c("CD14", "S100A8", "FCGR3A"),
  NK.cell = c("NCAM1", "NKG7", "TRDC"),
  B.cell = c("CD79A", "CD83", "MS4A1"),
  Classical.Dendritic = c("CLEC9A", "FCER1A", "CLEC10A"),
  Plasmacytoid.Dendritic = c("CLEC4C", "IL3RA", "LILRA4"),
  Plasma.cell = c("JCHAIN", "TNFRSF17", "SDC1")
)

# Visualize marker expression
for (g in names(marker.genes)) {
  FeaturePlot(seurat, features = marker.genes[[g]], order = TRUE, ncol = 3)
}

# Annotate clusters with cell type labels
seurat$celltype <- as.character(seurat$seurat_clusters)
seurat$celltype[seurat$celltype %in% c("0")] <- "NK.cell"
seurat$celltype[seurat$celltype %in% c("1", "5", "7", "8")] <- "T.cell"
seurat$celltype[seurat$celltype %in% c("2", "6")] <- "B.cell"
seurat$celltype[seurat$celltype %in% c("13")] <- "Plasma.cell"
seurat$celltype[seurat$celltype %in% c("3", "4", "9", "11")] <- "Monocyte"
seurat$celltype[seurat$celltype %in% c("10")] <- "Classical.Dendritic"
seurat$celltype[seurat$celltype %in% c("12")] <- "Plasmacytoid.Dendritic"

# Plot the annotated clusters
DimPlot(seurat, group.by = "celltype", label = FALSE)
