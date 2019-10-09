library(Seurat)
library(R.utils)
library(parallel)
library(future)

seed <- 0
n.cores <- detectCores()
logfile <- "seurat.log"
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores for tSNE."))

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan("multiprocess", workers = n.cores)
plan()

#### Benchmark Highly Variable Gene selection. ##
##src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_filter_norm.h5ad")
##obj.list <- SplitObject(src.obj, split.by = "Channel")
##print("Finding Highly Variable Genes:")
##now <- Sys.time()
##features <- SelectIntegrationFeatures(obj.list, nfeatures = 2000, fvf.nfeatures = 2000, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
##logstr.hvgsif <- Sys.time() - now
##write(paste("HVG + SIF time:", logstr.hvgsif, attr(logstr.hvgsif, "units")), file = logfile)
##
##src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_filter_norm.h5ad")
##obj.list <- SplitObject(src.obj, split.by = "Channel")
##print("HVG only:")
##now <- Sys.time()
##obj.list <- lapply(X = obj.list, FUN = function(x) {
##	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
##})
##logstr.hvg <- Sys.time() - now
##print(logstr.hvg)
##write(paste("HVG only time:", logstr.hvg, attr(logstr.hvg, "units")), file = logfile, append = TRUE)
##
#### Benchmark Batch Correction. ##
####src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_filter_norm_hvf.h5ad")
####VariableFeatures(src.obj) <- GetAssayData(src.obj)@Dimnames[[1]]
####obj.list <- SplitObject(src.obj, split.by = "Channel")
####now <- Sys.time()
####bm.anchors <- FindIntegrationAnchors(obj.list)
####bm.combined <- IntegrateData(bm.anchors)
####print("Batch Correction:")
####print(Sys.time() - now)
##
##
#### Benchmark PCA.
##adata <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_corrected_hvf.h5ad")
##print("PCA time:")
##now <- Sys.time()
##adata <- ScaleData(adata, scale.max = 10, verbose = FALSE)
##adata <- RunPCA(adata, verbose = FALSE, seed.use = seed, features = GetAssayData(adata, slot = "data")@Dimnames[[1]])
##logstr.pca <- Sys.time() - now
##print(logstr.pca)
##write(paste("PCA time:", logstr.pca, attr(logstr.pca, "units")), file = logfile, append = TRUE)
##rm(adata)

##load("MantonBM_nonmix.RData")
##print("Computing neighborhood graph:")
##n.pc <- dim(adata[["pca"]])[2]
##write(paste0("KNN using ", n.pc, " PCs"), file = logfile, append = TRUE)
##now <- Sys.time()
##adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy")
##logstr.knn <- Sys.time() - now
##write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile, append = TRUE)
##rm(adata)

##print("Finding Clusters using Louvain:")
##now <- Sys.time()
##adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "louvain")
##logstr.louvain <- Sys.time() - now
##write(paste("Louvain time:", logstr.louvain, attr(logstr.louvain, "units")), file = logfile, append = TRUE)

load("seurat_knn.RData")
library(Seurat)
source("/opt/software/seurat-3.1.0/R/clustering.R")
library(leiden)
library(igraph)

plan("multiprocess", workers = n.cores)
plan()

print("Finding Clusters using Leiden:")
now <- Sys.time()
#clustering.results <- FindClusters(adata[["RNA_snn"]], resolution = 1.3, algorithm = "leiden", method = "igraph", random.seed = seed)
input <- graph_from_adjacency_matrix(adjmatrix = adata[["RNA_snn"]])
ids <- leiden(
	object = input, 
	partition_type = "RBConfigurationVertexPartition",
    initial_membership = NULL,
    weights = NULL,
    node_sizes = NULL,
    resolution_parameter = 1.3,
    seed = seed,
    n_iterations = 10
)
names(x = ids) <- colnames(x = adata[["RNA_snn"]])
ids <- GroupSingletons(ids = ids, SNN = adata[["RNA_snn"]], group.singletons = TRUE, verbose = TRUE)
clustering.results[, paste0("res.", 1.3)] <- factor(x = ids)
colnames(x = clustering.results) <- paste0("RNA_snn", "_", colnames(x = clustering.results))
adata <- AddMetaData(object = adata, metadata = clustering.results)
Idents(object = adata) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
levels <- levels(x = adata)
Idents(object = adata) <- factor(x = Idents(object = adata), levels = sort(x = levels))
adata[['seurat_clusters']] <- Idents(object = adata)
adata <- LogSeuratCommand(adata)
logstr.leiden <- Sys.time() - now
write(paste("Leiden time:", logstr.leiden, attr(logstr.leiden, "units")), file = logfile, append = TRUE)

save(adata, file = "seurat_leiden.RData")

##
##load("MantonBM_nonmix.RData")
##n.pc <- dim(adata[["pca"]])[2]
##print("Computing UMAP embedding:")
##now <- Sys.time()
##adata <- RunUMAP(adata, reduction = "pca", umap.method = "umap-learn", n.neighbors = 15, min.dist = 0.5, seed.use = seed, dims = 1:n.pc)
##logstr.umap <- Sys.time() - now
##write(paste("UMAP time:", logstr.umap, attr(logstr.umap, "units")), file = logfile, append = TRUE)
##
##print("Computing FIt-SNE embedding:")
##now <- Sys.time()
##adata <- RunTSNE(adata, reduction = "pca", seed.use = seed, tsne.method = "FIt-SNE", dims = 1:n.pc)
##logstr.tsne <- Sys.time() - now
##print(logstr.tsne)
##write(paste("TSNE time:", logstr.tsne, attr(logstr.tsne, "units")), file = logfile, append = TRUE)

##save(adata, file = "seurat_result.RData")
