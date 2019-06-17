library(Seurat)
library(R.utils)
library(Matrix)
library(hash)
library(stringi)
library(ggplot2)
library(RANN)

seed <- 0
setOption('mc.cores', 8)

processed <- FALSE

if (!processed) {

    adata <- Read10X_h5("../MantonBM_nonmix_10x.h5")

    n.tags <- adata@Dim[2]
    subset.idx.hashmap <- hash()

    ##Split raw matrix into matrices due to channel.
    now <- Sys.time()
    for (j in 1:n.tags) {
        name <- adata@Dimnames[[2]][j]
        if (stri_detect(name, regex = "_HiSeq_1")) {
            chan <- strsplit(name, '-')[[1]][1]
            if (has.key(chan, subset.idx.hashmap)) {
                subset.idx.hashmap[[chan]] <- append(subset.idx.hashmap[[chan]], j)
            } else {
                subset.idx.hashmap[[chan]] <- c(j)
            }
        }
    }

## Split raw matrix into matrices due to individual.
##for (j in 1:n.tags) {
##    name <- adata@Dimnames[[2]][j]
##    person <- strsplit(name, '_')[[1]][1]
##    if (has.key(person, subset.idx.hashmap)) {
##        subset.idx.hashmap[[person]] <- append(subset.idx.hashmap[[person]], j)
##    } else {
##        subset.idx.hashmap[[person]] <- c(j)
##    }
##}

    dsets <- list()
    for (k in keys(subset.idx.hashmap)) {
        dsets[[k]] <- adata[, subset.idx.hashmap[[k]]]
    }
    print("Splitting 10X data:")
    print(Sys.time() - now)

## Read preset high variable features
##features <- as.character(unlist(read.table("hvg.txt")))
    df.hvg <- read.csv(file = "../hvg.txt")
    hvg.features <- as.character(unlist(df.hvg['index']))

## Create Seurat objects
    now <- Sys.time()
    seurat.obj.list <- list()
    for (i in 1:length(dsets)) {
        name <- names(dsets)[i]
        print(paste(i, ". Processing", name))
        ica <- CreateSeuratObject(counts = dsets[[i]], project = name, min.cells = 138)
        mito.features <- grep(pattern = "^MT-", rownames(ica), value = TRUE)
        percent.mito <- Matrix::colSums(GetAssayData(ica, slot = 'counts')[mito.features,]) / Matrix::colSums(GetAssayData(ica, slot = 'counts'))
        ica <- subset(ica, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mito < 0.1)
        ica <- NormalizeData(ica, selection.method = "LogNormalize", scale.factor = 1e5)
        #ica <- FindVariableFeatures(ica, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
        ica[["RNA"]]@var.features <- hvg.features
        seurat.obj.list[[name]] <- ica
    }
    print("Creating Seurat Objects:")
    print(Sys.time() - now)



    now <- Sys.time()
    ica.anchors <- FindIntegrationAnchors(object.list = seurat.obj.list)
    print("Computing Anchors using CCA:")
    print(Sys.time() - now)

    now <- Sys.time()
    ## Yiming: bug -- some features cannot find their indices.
    ica.combined <- IntegrateData(anchorset = ica.anchors)
    print("Data Integration:")
    print(Sys.time() - now)

    DefaultAssay(object = ica.combined) <- "integrated"

##print("Writing results to file...")
##myobj <- ica.combined[["integrated"]]
##X <- GetAssayData(object = myobj, slot = "data")
##writeMM(X, file = "count_matrix.mtx")
##write.table(X@Dimnames[[1]], file = "features.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
##write.table(X@Dimnames[[2]], file = "cells.txt", quote = FALSE, row.names= FALSE, col.names = FALSE)
##print("Finished!")

    save(ica.combined, file = "seurat_cca_processed.RData")
    print("Saved CCA result to file.")
} else {
    load(file = "seurat_cca_result.RData")
    print("Pre-calculated CCA result is loaded.")
}

now <- Sys.time()
ica.combined <- ScaleData(ica.combined, features = VariableFeatures(ica.combined), scale.max = 10)
print("Scale expression matrix of highly variable genes")

now <- Sys.time()
ica.combined <- RunPCA(ica.combined, features = VariableFeatures(ica.combined))
print("Running PCA:")
print(Sys.time() - now)

now <- Sys.time()
ica.combined <- FindNeighbors(ica.combined, k.param = 100, dims = 1:ncol(ica.combined[["pca"]]))
print("Computing neighborhood graph:")
print(Sys.time() - now)

now <- Sys.time()
ica.combined <- FindClusters(ica.combined, resolution = 1.3)
print("Finding Clusters:")
print(Sys.time() - now)

now <- Sys.time()
ica.combined <- RunUMAP(ica.combined, dims = 1:2, random_state = seed)
print("Computing UMAP embedding:")
print(Sys.time() - now)

now <- Sys.time()
ica.combined <- RunTSNE(ica.combined, seed.use = seed, tsne.method = "FIt-SNE")
print("Computing tSNE embedding:")
print(Sys.time() - now)

save(ica.combined, file = "seurat_final.RData")
print("Saved final result to file.")

## Plot TSNE and UMAP
p1 <- DimPlot(ica.combined, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
p2 <- DimPlot(ica.combined, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
CombinePlots(plots = list(p1, p2), legend = "none")
ggsave(filename = "seurat_cca_fitsne+umap.pdf", width = 20, height = 10)

## Save PCA and Cell information
pca.obj <- ica.combined[["pca"]]
X.pca <- pca.obj@cell.embeddings
write.table(X.pca, file = "seurat_cca_pca.txt", quote = FALSE, sep = ",", row.names = FALSE)

cells <- unlist(dimnames(pca.obj@cell.embeddings)[[1]])
write.table(cells, file = "seurat_cca_cells.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

## Calculate KNN indices, and save to file.
my.knn <- nn2(data = X.pca, k = 100, searchtype = "standard", eps = 0)
write.table(my.knn$nn.idx, "seurat_cca_knn.txt", row.names = FALSE, sep = ",", quote = FALSE)
