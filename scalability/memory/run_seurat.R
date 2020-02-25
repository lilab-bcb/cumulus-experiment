library(Seurat)
library(R.utils)
library(parallel)
library(future)

record.time <- function(start.time, end.time, title, filename) {
	duration.time <- end.time - start.time
    write(paste(title, "starts at:", start.time), file = filename, append = TRUE)
    write(paste(title, "ends at:", end.time), file = filename, append = TRUE)
    write(paste(title, "time:", duration.time, attr(duration.time, "units")), file = filename, append = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)

src.data <- args[2]
print(src.data)
out.name <- args[3]
print(out.name)

seed <- 0
n.cores <- detectCores()
logfile <- paste0(out.name, "_seurat.log")
debug.mode <- FALSE
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores for tSNE."))
write(paste("Use", n.cores, "cores for tSNE."), file = logfile)

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan("multiprocess", workers = n.cores)
plan()

start.read <- Sys.time()
src.obj <- Read10X_h5(src.data)
end.read <- Sys.time()
record.time(start.read, end.read, "Read", logfile)

start.filter <- Sys.time()
n.cells <- dim(src.obj)[2]
bm <- CreateSeuratObject(src.obj, project = "MantonBM", min.cells = floor(n.cells * 0.0005), min.features = 500)
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm <- subset(bm, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mt < 10)
end.filter <- Sys.time()
record.time(start.filter, end.filter, "Filter", logfile)

start.norm <- Sys.time()
bm <- NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 1e5, verbose = debug.mode)
end.norm <- Sys.time()
record.time(start.norm, end.norm, "Normalization", logfile)

start.hvg <- Sys.time()
bm <- FindVariableFeatures(bm, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000, verbose = debug.mode)
end.hvg <- Sys.time()
record.time(start.hvg, end.hvg, "HVG", logfile)

start.pca <- Sys.time()
bm <- ScaleData(bm, scale.max = 10, verbose = debug.mode)
bm <- RunPCA(bm, npcs = 50, seed.use = seed, verbose = debug.mode)
end.pca <- Sys.time()
record.time(start.pca, end.pca, "PCA", logfile)

start.knn <- Sys.time()
bm <- FindNeighbors(bm, reduction = "pca", dims = 1:50, nn.method = 'annoy', k.param = 100, verbose = debug.mode)
end.knn <- Sys.time()
record.time(start.knn, end.knn, "KNN", logfile)

start.louvain <- Sys.time()
bm <- FindClusters(bm, resolution = 1.3, algorithm = "louvain", random.seed = seed, verbose = debug.mode)
end.louvain <- Sys.time()
record.time(start.louvain, end.louvain, "Louvain", logfile)

start.tsne <- Sys.time()
bm <- RunTSNE(bm, reduction = "pca", dims = 1:50, tsne.method = 'FIt-SNE', seed.use = seed)
end.tsne <- Sys.time()
record.time(start.tsne, end.tsne, "tSNE", logfile)

start.umap <- Sys.time()
bm <- RunUMAP(bm, reduction = "pca", umap.method = "umap-learn", metric = "correlation", n.neighbors = 15, dims = 1:50, seed.use = seed, verbose = debug.mode)
end.umap <- Sys.time()
record.time(start.umap, end.umap, "UMAP", logfile)

save(bm, file = paste0(out.name, "_seurat_result.RData"))