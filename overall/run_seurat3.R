library(Seurat)
library(R.utils)
library(parallel)

seed <- 0
n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))

## Benchmark Highly Variable Gene selection. ##
##src.obj <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm.h5ad")
##obj.list <- SplitObject(src.obj, split.by = "Channel")
##now <- Sys.time()
##hvg <- SelectIntegrationFeatures(obj.list, nfeatures = 2000, fvf.nfeatures = 2000, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
##print("Finding Highly Variable Genes:")
##print(Sys.time() - now)

## Benchmark Batch Correction. ##
src.obj <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm_hvg.h5ad")
VariableFeatures(src.obj) <- GetAssayData(src.obj)@Dimnames[[1]]
obj.list <- SplitObject(src.obj, split.by = "Channel")
now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(obj.list)
bm.combined <- IntegrateData(bm.anchors)
print("Batch Correction:")
print(Sys.time() - now)



if (1 == 0) {

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

print("Finding Clusters:")
now <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain")
print("Clustering time:")
print(Sys.time() - now)

save(adata, file = "seurat_cluster.RData")

print("Computing UMAP embedding:")
now <- Sys.time()
ica <- RunUMAP(ica, n.neighbors = 15, min.dist = 0.5, seed.use = seed, features = ica[["pca"]], metric = 'euclidean')
print("UMAP time:")
print(Sys.time() - now)

print("Computing FItSNE embedding:")
now <- Sys.time()
ica <- RunTSNE(ica, seed.use = seed, tsne.method = "FIt-SNE", features = ica[["pca"]])
print("FItSNE time:")
print(Sys.time() - now)

}
