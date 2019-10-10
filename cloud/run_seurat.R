library(Seurat)
library(R.utils)
library(parallel)
library(future)

src_file <- "/projects/benchmark/MantonBM/MantonBM_nonmix_10x.h5"
dst_folder <- "./"
seed <- 0
debug.mode <- TRUE
logfile <- paste0(dst_folder, "seurat.log")

n.cores <- detectCores()
setOption('mc.cores', n.cores)
write(paste("Use", n.cores, "cores for tSNE."), file = logfile)

options(future.globals.maxSize = 20 * 1024^3)
plan("multiprocess", workers = n.cores)
plan()
write(paste("Use", nbrOfWorkers(), "cores."), file = logfile, append = TRUE)

bm.data <- Read10X_h5(src_file)
n.cells <- dim(bm.data)[2]
bm <- CreateSeuratObject(bm.data, project = "MantonBM", min.cells = floor(n.cells * 0.0005), min.features = 500)
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm <- subset(bm, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mt < 10)

bm.list <- SplitObject(bm)

bm.list <- lapply(X = bm.list, FUN = function(x) {
	x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e5, verbose = debug.mode)
})

now <- Sys.time()
bm.list <- lapply(X = bm.list, FUN = function(x) {
	x <- FindVariableFeatures(x, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000, verbose = debug.mode)
})
logstr.hvg <- Sys.time() - now
write(paste("HVG time:", logstr.hvg, attr(logstr.hvg, "units")), file = logfile, append = TRUE)

plan("multiprocess", workers = 15)
plan()
write(paste("Use", nbrOfWorkers(), "cores."), file = logfile, append = TRUE)

now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, dims = 1:20, verbose = debug.mode)
logstr.anchor <- Sys.time() - now
write(paste("Find Anchors time:", logstr.anchor, attr(logstr.anchor, "units")), file = logfile, append = TRUE)

save(bm.anchors, file = paste0(dst_folder, "seurat_anchor_result.RData"))
rm(bm.list)

now <- Sys.time()
bm.combined <- IntegrateData(anchorset = bm.anchors, dims = 1:20, verbose = debug.mode)
logstr.cca <- Sys.time() - now
write(paste("Integration time:", logstr.cca, attr(logstr.cca, "units")), file = logfile, append = TRUE)

DefaultAssay(bm.combined) <- "integrated"

save(bm.combined, file = paste0(dst_folder, "seurat_corrected_result.RData"))
rm(bm.anchors)

plan("multiprocess", workers = n.cores)
plan()
write(paste("Use", nbrOfWorkers(), "cores."), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- ScaleData(bm.combined, scale.max = 10, verbose = debug.mode)
bm.combined <- RunPCA(bm.combined, npcs = 50, seed.use = seed, verbose = debug.mode)
logstr.pca <- Sys.time() - now
write(paste("PCA time:", logstr.pca, attr(logstr.pca, "units")), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- FindNeighbors(bm.combined, reduction = "pca", dims = 1:50, nn.method = 'annoy', k.param = 100, compute.SNN = TRUE, verbose = debug.mode)
logstr.knn <- Sys.time() - now
write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- FindClusters(bm.combined, resolution = 1.3, algorithm = "louvain", random.seed = seed, verbose = debug.mode)
logstr.louvain <- Sys.time() - now
write(paste("Louvain time:", logstr.louvain, attr(logstr.louvain, "units")), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- RunTSNE(bm.combined, reduction = "pca", dims = 1:50, tsne.method = 'FIt-SNE', seed.use = seed)
logstr.tsne <- Sys.time() - now
write(paste("TSNE time:", logstr.tsne, attr(logstr.tsne, "units")), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- RunUMAP(bm.combined, reduction = "pca", umap.method = "umap-learn", metric = "correlation", n.neighbors = 15, dims = 1:50, seed.use = seed, verbose = debug.mode)
logstr.umap <- Sys.time() - now
write(paste("UMAP time:", logstr.umap, attr(logstr.umap, "units")), file = logfile, append = TRUE)

DefaultAssay(bm.combined) <- "RNA"
save(bm.combined, file = paste0(dst_folder, "seurat_pipeline_result.RData"))

now <- Sys.time()
bm.markers <- FindAllMarkers(bm.combined, test.use = "t", random.seed = seed, verbose = debug.mode)
logstr.de <- Sys.time() - now
write(paste("DE time:", logstr.de, attr(logstr.de, "units")), file = logfile, append = TRUE)

save(bm.markers, file = paste0(dst_folder, "seurat_de_analysis_result.RData"))