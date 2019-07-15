library(Seurat)
library(R.utils)
library(parallel)

seed <- 0
n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))


## Benchmark Highly Variable Gene selection.
ica <- ReadH5AD("../MantonBM_nonmix_10x_filter_norm.h5ad")
now <- Sys.time()
ica <- FindVariableFeatures(ica, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
print("Finding Highly Variable Genes:")
print(Sys.time() - now)

## Benchmark PCA.
##adata <- ReadH5AD("./sccloud_output/MantonBM_nonmix_sccloud_corrected.h5ad")
##now <- Sys.time()
####ica <- ScaleData(ica, features = rownames(ica), scale.max = 10)
##ica <- ScaleData(ica, scale.max = 10)
##print("Scale expression matrix of variable genes:")
##print(Sys.time() - now)
##
##now <- Sys.time()
##print("Running PCA:")
##ica <- RunPCA(ica, features = VariableFeatures(ica), verbose = FALSE, seed.use = seed)
##print("PCA time:")
##print(Sys.time() - now)
##
##now <- Sys.time()
##print("Computing neighborhood graph:")
##ica <- FindNeighbors(ica, k.param = 100, compute.SNN = FALSE) # Issue: Do I have to compute SNN?
##print("kNN time:")
##print(Sys.time() - now)
##
##now <- Sys.time()
##print("Finding Clusters using Leiden:")
##ica <- FindClusters(ica, resolution = 1.3, random.seed = seed, algorithm = 4, verbose = FALSE) # TODO: Check if 2 is faster than 4.
##print("Leiden Clustering time:")
##print(Sys.time() - now)
##
##now <- Sys.time()
##print("Computing UMAP embedding:")
##ica <- RunUMAP(ica, min.dist = 0.5, seed.use = seed, features = ica[["pca"]], verbose = FALSE)
##print("UMAP time:")
##print(Sys.time() - now)
##
##now <- Sys.time()
##print("Computing FItSNE embedding:")
##ica <- RunTSNE(ica, seed.use = seed, tsne.method = "FIt-SNE", features = ica[["pca"]])
##print("FItSNE time:")
##print(Sys.time() - now)