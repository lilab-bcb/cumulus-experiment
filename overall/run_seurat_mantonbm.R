library(Seurat)
library(R.utils)
library(parallel)
library(future)

seed <- 0
n.cores <- detectCores()
logfile <- "mantonbm_seurat.log"
debug.mode <- FALSE
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores for tSNE."))

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan("multiprocess", workers = n.cores)
plan()

## Benchmark Highly Variable Gene selection. ##
src.obj <- ReadH5AD("/data/MantonBM_nonmix_filter_norm.h5ad")
obj.list <- SplitObject(src.obj, split.by = "Channel")
print("Finding Highly Variable Genes:")
start.hvg <- Sys.time()
features <- SelectIntegrationFeatures(obj.list, verbose = debug.mode, nfeatures = 2000, fvf.nfeatures = 2000, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
end.hvg <- Sys.time()
duration.hvg <- end.hvg - start.hvg
write(paste("HVG starts at:", start.hvg), file = logfile)
write(paste("HVG ends at:", end.hvg), file = logfile, append = TRUE)
write(paste("HVG time:", duration.hvg, attr(duration.hvg, "units")), file = logfile, append = TRUE)

## Benchmark PCA.
adata <- ReadH5AD("/data/MantonBM_nonmix_corrected_hvf.h5ad")
print("PCA:")
start.pca <- Sys.time()
adata <- ScaleData(adata, scale.max = 10, verbose = debug.mode)
adata <- RunPCA(adata, verbose = debug.mode, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
end.pca <- Sys.time()
duration.pca <- end.pca - start.pca
write(paste("PCA starts at:", start.pca), file = logfile, append = TRUE)
write(paste("PCA ends at:", end.pca), file = logfile, append = TRUE)
write(paste("PCA time:", duration.pca, attr(duration.pca, "units")), file = logfile, append = TRUE)
rm(adata)

## Benchmark KNN.
load("MantonBM_nonmix.RData")
print("Computing neighborhood graph:")
n.pc <- dim(adata[["pca"]])[2]
start.knn <- Sys.time()
bdata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = FALSE, verbose = debug.mode)
end.knn <- Sys.time()
duration.knn <- end.knn - start.knn
write(paste("KNN starts at:", start.knn), file = logfile, append = TRUE)
write(paste("KNN ends at:", end.knn), file = logfile, append = TRUE)
write(paste("KNN time:", duration.knn, attr(duration.knn, "units")), file = logfile, append = TRUE)
rm(bdata)

adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = TRUE, verbose = debug.mode)

print("Finding Clusters using Louvain:")
start.louvain <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain", verbose = debug.mode)
end.louvain <- Sys.time()
duration.louvain <- end.louvain - start.louvain
write(paste("Louvain starts at:", start.louvain), file = logfile, append = TRUE)
write(paste("Louvain ends at:", end.louvain), file = logfile, append = TRUE)
write(paste("Louvain time:", duration.louvain, attr(duration.louvain, "units")), file = logfile, append = TRUE)
rm(adata)

load("MantonBM_nonmix.RData")
n.pc <- dim(adata[["pca"]])[2]
print("Computing UMAP embedding:")
start.umap <- Sys.time()
adata <- RunUMAP(adata, reduction = "pca", umap.method = "umap-learn", metric = "correlation", n.neighbors = 15, min.dist = 0.5, seed.use = seed, dims = 1:n.pc, verbose = debug.mode)
end.umap <- Sys.time()
duration.umap <- end.umap - start.umap
write(paste("UMAP starts at:", start.umap), file = logfile, append = TRUE)
write(paste("UMAP ends at:", end.umap), file = logfile, append = TRUE)
write(paste("UMAP time:", duration.umap, attr(duration.umap, "units")), file = logfile, append = TRUE)

print("Computing FIt-SNE embedding:")
start.tsne <- Sys.time()
adata <- RunTSNE(adata, reduction = "pca", seed.use = seed, tsne.method = "FIt-SNE", dims = 1:n.pc)
end.tsne <- Sys.time()
duration.tsne <- end.tsne - start.tsne
write(paste("TSNE starts at:", start.tsne), file = logfile, append = TRUE)
write(paste("TSNE ends at:", end.tsne), file = logfile, append = TRUE)
write(paste("TSNE time:", duration.tsne, attr(duration.tsne, "units")), file = logfile, append = TRUE)