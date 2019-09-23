library(Seurat)
library(R.utils)
library(parallel)

seed <- 0
n.cores <- detectCores()
logfile <- "seurat.log"
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))

## Benchmark Highly Variable Gene selection. ##
src.obj <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm.h5ad")
obj.list <- SplitObject(src.obj, split.by = "Channel")
print("Finding Highly Variable Genes:")
now <- Sys.time()
hvg <- SelectIntegrationFeatures(obj.list, nfeatures = 2000, fvf.nfeatures = 2000, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
logstr.hvg <- Sys.time() - now
print(logstr.hvg)
write(paste("HVG time:", logstr.hvg, attr(logstr.hvg, "units")), file = logfile)

## Benchmark Batch Correction. ##
##src.obj <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm_hvg.h5ad")
##VariableFeatures(src.obj) <- GetAssayData(src.obj)@Dimnames[[1]]
##obj.list <- SplitObject(src.obj, split.by = "Channel")
##now <- Sys.time()
##bm.anchors <- FindIntegrationAnchors(obj.list)
##bm.combined <- IntegrateData(bm.anchors)
##print("Batch Correction:")
##print(Sys.time() - now)


## Benchmark PCA.
adata <- ReadH5AD("../MantonBM_nonmix_corrected_hvg.h5ad")
print("PCA time:")
now <- Sys.time()
adata <- ScaleData(adata, scale.max = 10, verbose = FALSE)
adata <- RunPCA(adata, verbose = FALSE, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
logstr.pca <- Sys.time() - now
print(logstr.pca)
write(paste("PCA time:", logstr.pca, attr(logstr.pca, "units")), file = logfile, append = TRUE)


save(adata, file = "seurat_pca.RData")

print("Computing neighborhood graph:")
n.pc <- dim(adata[["pca"]])[2]
print(paste0("KNN using ", n.pc, " PCs."))
write(paste0("KNN using ", n.pc, " PCs"), file = logfile, append = TRUE)
now <- Sys.time()
bdata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = FALSE)
logstr.knn <- Sys.time() - now
print(logstr.knn)
write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile, append = TRUE)

adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = TRUE)
save(adata, file = "seurat_knn.RData")


print("Finding Clusters:")
now <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain")
logstr.cluster <- Sys.time() - now
print(logstr.cluster)
write(paste("Clustering time:", logstr.cluster, attr(logstr.cluster, "units")), file = logfile, append = TRUE)

print("Computing UMAP embedding:")
now <- Sys.time()
adata <- RunUMAP(adata, n.neighbors = 15, min.dist = 0.5, seed.use = seed, dims = 1:n.pc)
logstr.umap <- Sys.time() - now
print(logstr.umap)
write(paste("UMAP time:", logstr.umap, attr(logstr.umap, "units")), file = logfile, append = TRUE)


print("Computing FItSNE embedding:")
now <- Sys.time()
adata <- RunTSNE(adata, seed.use = seed, tsne.method = "FIt-SNE", dims = 1:n.pc)
logstr.tsne <- Sys.time() - now
print(logstr.tsne)
write(paste("TSNE time:", logstr.tsne, attr(logstr.tsne, "units")), file = logfile, append = TRUE)
