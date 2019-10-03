library(Seurat)
library(R.utils)
library(parallel)
library(future)

seed <- 0
n.cores <- detectCores()
logfile <- "seurat.log"
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores for tSNE."))

plan("multiprocess", workers = 32)
plan()

bm.data <- Read10X_h5("/projects/benchmark/MantonBM/MantonBM_nonmix_10x.h5")
n.cells <- dim(bm.data)[2]
bm <- CreateSeuratObject(bm.data, project = "MantonBM", min.cells = floor(n.cells * 0.0005), min.features = 500)
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm <- subset(bm, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mt < 10)


bm.list <- SplitObject(bm)
bm.list <- lapply(X = bm.list, FUN = function(x) {
	x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e5)
	now <- Sys.time()
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
	logstr.hvg <- Sys.time() - now
	write(paste("HVG time:", logstr.hvg, attr(logstr.hvg, "units")), file = logfile)
})

now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, dims = 1:20)
bm.combined <- IntegrateData(anchorset = bm.anchors, dims = 1:20)
logstr.cca <- Sys.time() - now
write(paste("Batch Correction time:", logstr.cca, attr(logstr.cca, "units")), file = logfile)

DefaultAssay(bm.combined) <- "integrated"

now <- Sys.time()
bm.combined <- ScaleData(bm.combined)
bm.combined <- RunPCA(bm.combined, npcs = 50, seed.use = seed)
logstr.pca <- Sys.time() - now
write(paste("PCA time:", logstr.pca, attr(logstr.pca, "units")), file = logfile)

now <- Sys.time()
bm.combined <- RunTSNE(bm.combined, reduction = "pca", dims = 1:50, tsne.method = 'FIt-SNE', seed.use = seed)
logstr.tsne <- Sys.time() - now
write(paste("TSNE time:", logstr.tsne, attr(logstr.tsne, "units")), file = logfile)

now <- Sys.time()
bm.combined <- RunUMAP(bm.combined, reduction = "pca", dims = 1:50, seed.use = seed)
logstr.umap <- Sys.time() - now
write(paste("UMAP time:", logstr.umap, attr(logstr.umap, "units")), file = logfile)

now <- Sys.time()
bm.combined <- FindNeighbors(bm.combined, reduction = "pca", dims = 1:50, nn.method = 'annoy', k.param = 100)
logstr.knn <- Sys.time() - now
write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile)

now <- Sys.time()
bm.combined <- FindClusters(bm.combined, resolution = 1.3, algorithm = 1, random.seed = 0)
logstr.louvain <- Sys.time() - now
write(paste("Louvain time:", logstr.louvain, attr(logstr.louvain, "units")), file = logfile)

DefaultAssay(bm.combined) <- "RNA"
now <- Sys.time()
bm.markers <- FindAllMarkers(bm.combined, test.use = "t", random.seed = 0)
logstr.de <- Sys.time() - now
write(paste("DE time:", logstr.de, attr(logstr.de, "units")), file = logfile)