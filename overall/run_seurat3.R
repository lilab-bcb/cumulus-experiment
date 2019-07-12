library(Seurat)
library(R.utils)
library(parallel)

seed <- 0
n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))

start <- Sys.time()
adata <- Read10X_h5("../MantonBM_nonmix_10x.h5")
rownames(adata) <- make.names(adata[,1], unique = TRUE)
ica <- CreateSeuratObject(counts = adata, project = "ICA", min.cells = 138)
print("Read and Initialization:")

now <- Sys.time()
mito.features <- grep(pattern = "^MT-", rownames(ica), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(ica, slot = 'counts')[mito.features,]) / Matrix::colSums(GetAssayData(ica, slot = 'counts'))
ica <- subset(ica, subset = nFeature_RNA >=500 & nFeature_RNA < 6000 & percent.mito < 0.1)
print("Filtering cells:")
print(Sys.time() - now)

now <- Sys.time()
ica <- NormalizeData(ica, selection.method = "LogNormalize", scale.factor = 1e5)
print("Normalizing data:")
print(Sys.time() - now)

now <- Sys.time()
ica <- FindVariableFeatures(ica, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
print("Finding variable genes:")
print(Sys.time() - now)

now <- Sys.time()
##ica <- ScaleData(ica, features = rownames(ica), scale.max = 10)
ica <- ScaleData(ica, scale.max = 10)
print("Scale expression matrix of variable genes:")
print(Sys.time() - now)

now <- Sys.time()
print("Running PCA:")
ica <- RunPCA(ica, features = VariableFeatures(ica), verbose = FALSE, seed.use = seed)
print("PCA time:")
print(Sys.time() - now)

now <- Sys.time()
print("Computing neighborhood graph:")
ica <- FindNeighbors(ica, k.param = 100, compute.SNN = FALSE) # Issue: Do I have to compute SNN?
print("kNN time:")
print(Sys.time() - now)

now <- Sys.time()
print("Finding Clusters using Leiden:")
ica <- FindClusters(ica, resolution = 1.3, random.seed = seed, algorithm = 4, verbose = FALSE) # TODO: Check if 2 is faster than 4.
print("Leiden Clustering time:")
print(Sys.time() - now)

now <- Sys.time()
print("Computing UMAP embedding:")
ica <- RunUMAP(ica, min.dist = 0.5, seed.use = seed, features = ica[["pca"]], verbose = FALSE)
print("UMAP time:")
print(Sys.time() - now)

now <- Sys.time()
print("Computing FItSNE embedding:")
ica <- RunTSNE(ica, seed.use = seed, tsne.method = "FIt-SNE", features = ica[["pca"]])
print("FItSNE time:")
print(Sys.time() - now)

#now <- Sys.time()
#ica.markers <- FindAllMarkers(ica, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#print("Finding Markers:")
#print(Sys.time() - now)

end <- Sys.time()
print("Total time without writing:")
print(end - start)

#write.csv(ica.markers, file = "/home/de_analysis_seurat3.csv")
save(ica, file = "./seurat_output/seurat_result.RData")
end <- Sys.time()
print("Total time after R writing:")
print(end - start)