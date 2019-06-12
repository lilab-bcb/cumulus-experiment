library(Seurat)
library(R.utils)
library(Matrix)
library(hash)
library(stringi)
library(ggplot2)
library(RANN)

seed <- 0
setOption('mc.cores', 12)

adata <- Read10X_h5("../MantonBM_tiny_10x.h5")

## Create Seurat objects
ica <- CreateSeuratObject(counts = adata, project = "ICA", min.cells = 138)
mito.features <- grep(pattern = "^MT-", rownames(ica), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(ica, slot = 'counts')[mito.features,]) / Matrix::colSums(GetAssayData(ica, slot = 'counts'))
ica <- subset(ica, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mito < 0.1)

now <- Sys.time()
ica <- NormalizeData(ica, selection.method = "LogNormalize", scale.factor = 1e5)
print("Normalizing data:")
print(Sys.time() - now)

now <- Sys.time()
ica <- FindVariableFeatures(ica, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
print("Finding variable genes:")
print(Sys.time() - now)

now <- Sys.time()
ica <- ScaleData(ica, features = VariableFeatures(ica), scale.max = 10)
print("Scale expression matrix of highly variable genes:")
print(Sys.time() - now)

now <- Sys.time()
ica <- RunPCA(ica, features = VariableFeatures(ica), verbose = FALSE)
print("Running PCA:")
print(Sys.time() - now)

now <- Sys.time()
ica <- FindNeighbors(ica, k.param = 100, dims = 1:ncol(ica[["pca"]]))
print("Computing neighborhood graph:")
print(Sys.time() - now)

now <- Sys.time()
ica <- FindClusters(ica, resolution = 1.3)
print("Finding Clusters:")
print(Sys.time() - now)

now <- Sys.time()
ica <- RunUMAP(ica, dims = 1:2, random_state = seed)
print("Computing UMAP embedding:")
print(Sys.time() - now)

now <- Sys.time()
ica <- RunTSNE(ica, seed.use = seed, tsne.method = "FIt-SNE")
print("Computing tSNE embedding:")
print(Sys.time() - now)

save(ica, file = "seurat_raw_final.RData")
print("Saved final result to file.")

## Plot TSNE and UMAP
p1 <- DimPlot(ica, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
p2 <- DimPlot(ica, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
CombinePlots(plots = list(p1, p2), legend = "none")
ggsave(filename = "seurat_raw_fitsne+umap.pdf", width = 20, height = 10)

## Save PCA and Cell information
pca.obj <- ica[["pca"]]
X.pca <- pca.obj@cell.embeddings
write.table(X.pca, file = "seurat_raw_pca.txt", quote = FALSE, sep = ",", row.names = FALSE)

cells <- unlist(dimnames(pca.obj@cell.embeddings)[[1]])
write.table(cells, file = "seurat_raw_cells.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

## Calculate KNN indices, and save to file.
my.knn <- nn2(data = X.pca, k = 100, searchtype = "standard", eps = 0)
write.table(my.knn$nn.idx, file = "seurat_raw_knn.txt", row.names = FALSE, sep = ",", quote = FALSE)
