library(Seurat)
library(R.utils)
library(parallel)
library(future)

record.time <- function(duration.time, title, filename) {
    write(paste(title, "time:", duration.time, attr(duration.time, "units")), file = filename, append = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)

seed <- 0
n.cores <- detectCores()
if (length(args) == 2) {
    n.cores <- as.numeric(args[2])
}
logfile <- paste0("pbmc_seurat_cpu_", n.cores, ".log")
debug.mode <- FALSE
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores for tSNE."))
write(paste("Use", n.cores, "cores for tSNE."), file = logfile)

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan("multiprocess", workers = n.cores)
plan()

## Benchmark Highly Variable Gene selection. ##
src.obj <- ReadH5AD("/data/5k_pbmc_v3_filter_norm.h5ad")
print("Finding Highly Variable Genes:")
start.hvg <- Sys.time()
pbmc <- FindVariableFeatures(src.obj, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf), verbose = debug.mode)
end.hvg <- Sys.time()
record.time(end.hvg - start.hvg, "HVG", logfile)

## Benchmark PCA.
adata <- ReadH5AD("/data/5k_pbmc_v3_filter_norm_hvf.h5ad")
print("PCA:")
start.pca <- Sys.time()
adata <- ScaleData(adata, scale.max = 10, verbose = debug.mode)
adata <- RunPCA(adata, verbose = debug.mode, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
end.pca <- Sys.time()
record.time(end.pca - start.pca, "PCA", logfile)
rm(adata)

## Benchmark KNN.
load("5k_pbmc_v3.RData")
print("Computing neighborhood graph:")
n.pc <- dim(adata[["pca"]])[2]
start.knn <- Sys.time()
bdata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = FALSE, verbose = debug.mode)
end.knn <- Sys.time()
record.time(end.knn - start.knn, "KNN", logfile)
rm(bdata)

adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = TRUE, verbose = debug.mode)
save(adata, file = "pbmc_seurat_knn.RData")

print("Finding Clusters using Louvain:")
start.louvain <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain", verbose = debug.mode)
end.louvain <- Sys.time()
record.time(end.louvain - start.louvain, "Louvain", logfile)
save(adata, file = "pbmc_seurat_louvain.RData")
rm(adata)

load("pbmc_seurat_knn.RData")
print("Finding Clusters using Leiden:")
start.leiden <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "leiden", verbose = debug.mode)
end.leiden <- Sys.time()
record.time(end.leiden - start.leiden, "Leiden", logfile)
save(adata, file = "pbmc_seurat_leiden.RData")
rm(adata)

load("5k_pbmc_v3.RData")
n.pc <- dim(adata[["pca"]])[2]
print("Computing UMAP embedding:")
start.umap <- Sys.time()
adata <- RunUMAP(adata, reduction = "pca", umap.method = "umap-learn", metric = "correlation", n.neighbors = 15, min.dist = 0.5, seed.use = seed, dims = 1:n.pc, verbose = debug.mode)
end.umap <- Sys.time()
record.time(end.umap - start.umap, "UMAP", logfile)

print("Computing FIt-SNE embedding:")
start.tsne <- Sys.time()
adata <- RunTSNE(adata, reduction = "pca", seed.use = seed, tsne.method = "FIt-SNE", dims = 1:n.pc)
end.tsne <- Sys.time()
record.time(end.tsne - start.tsne, "tSNE", logfile)

save(adata, file = "pbmc_seurat_embedding.RData")
