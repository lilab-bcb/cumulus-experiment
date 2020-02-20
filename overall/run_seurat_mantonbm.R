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
now <- Sys.time()
features <- SelectIntegrationFeatures(obj.list, verbose = debug.mode, nfeatures = 2000, fvf.nfeatures = 2000, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
logstr.hvgsif <- Sys.time() - now
write(paste("HVG time:", logstr.hvgsif, attr(logstr.hvgsif, "units")), file = logfile)

## Benchmark PCA.
adata <- ReadH5AD("/data/MantonBM_nonmix_corrected_hvf.h5ad")
print("PCA:")
now <- Sys.time()
adata <- ScaleData(adata, scale.max = 10, verbose = debug.mode)
adata <- RunPCA(adata, verbose = debug.mode, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
logstr.pca <- Sys.time() - now
write(paste("PCA time:", logstr.pca, attr(logstr.pca, "units")), file = logfile, append = TRUE)
rm(adata)

## Benchmark KNN.
load("MantonBM_nonmix.RData")
print("Computing neighborhood graph:")
n.pc <- dim(adata[["pca"]])[2]
now <- Sys.time()
bdata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = FALSE, verbose = debug.mode)
logstr.knn <- Sys.time() - now
write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile, append = TRUE)

adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = TRUE, verbose = debug.mode)
save(adata, file = "mantonbm_seurat_knn.RData")

print("Finding Clusters using Louvain:")
now <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain", verbose = debug.mode)
logstr.louvain <- Sys.time() - now
write(paste("Louvain time:", logstr.louvain, attr(logstr.louvain, "units")), file = logfile, append = TRUE)
save(adata, file = "mantonbm_seurat_louvain.RData")
rm(adata)

load("MantonBM_nonmix.RData")
n.pc <- dim(adata[["pca"]])[2]
print("Computing UMAP embedding:")
now <- Sys.time()
adata <- RunUMAP(adata, reduction = "pca", umap.method = "umap-learn", metric = "correlation", n.neighbors = 15, min.dist = 0.5, seed.use = seed, dims = 1:n.pc, verbose = debug.mode)
logstr.umap <- Sys.time() - now
write(paste("UMAP time:", logstr.umap, attr(logstr.umap, "units")), file = logfile, append = TRUE)

print("Computing FIt-SNE embedding:")
now <- Sys.time()
adata <- RunTSNE(adata, reduction = "pca", seed.use = seed, tsne.method = "FIt-SNE", dims = 1:n.pc)
logstr.tsne <- Sys.time() - now
write(paste("TSNE time:", logstr.tsne, attr(logstr.tsne, "units")), file = logfile, append = TRUE)

save(adata, file = "mantonbm_seurat_visualization.RData")
