library(Seurat)
library(R.utils)
library(parallel)

seed <- 0
n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))

skip.finished <- TRUE

if (!skip.finished) {
## Benchmark Highly Variable Gene selection.
ica <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm.h5ad")
now <- Sys.time()
ica <- FindVariableFeatures(ica, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf), verbose = FALSE)
print("Finding Highly Variable Genes:")
print(Sys.time() - now)

}

## Benchmark PCA.
adata <- ReadH5AD("../MantonBM_nonmix_corrected_hvg.h5ad")
now <- Sys.time()
adata <- ScaleData(adata, scale.max = 10, verbose = FALSE)
adata <- RunPCA(adata, verbose = FALSE, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
print("PCA time:")
print(Sys.time() - now)

save(adata, file = "seurat_pca.RData")

print("Computing neighborhood graph:")
n.pc <- dim(adata[["pca"]])[2]
now <- Sys.time()
adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc)
print("kNN time in total:")
print(Sys.time() - now)

save(adata, file = "seurat_knn.RData")

print("Finding Clusters using Leiden:")
now <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = 4, verbose = FALSE)
print("Leiden Clustering time:")
print(Sys.time() - now)

save(adata, file = "seurat_cluster.RData")

print("Computing UMAP embedding:")
now <- Sys.time()
ica <- RunUMAP(ica, min.dist = 0.5, seed.use = seed, features = ica[["pca"]], verbose = FALSE)
print("UMAP time:")
print(Sys.time() - now)

print("Computing FItSNE embedding:")
now <- Sys.time()
ica <- RunTSNE(ica, seed.use = seed, tsne.method = "FIt-SNE", features = ica[["pca"]])
print("FItSNE time:")
print(Sys.time() - now)
