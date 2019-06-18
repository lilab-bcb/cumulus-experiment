library(Seurat)
library(R.utils)
library(Matrix)
library(hash)
library(stringi)
library(ggplot2)
library(RANN)

seed <- 0
setOption('mc.cores', 8)

src.obj <- ReadH5AD("../MantonBM_nonmix_tiny_10x_filter_norm.h5ad")
adata <- GetAssayData(src.obj)

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
df.hvg <- read.csv(file = "../hvg.txt")
hvg.features <- as.character(unlist(df.hvg['index']))

## Create Seurat objects
now <- Sys.time()
seurat.obj.list <- list()
for (i in 1:length(dsets)) {
    name <- names(dsets)[i]
    print(paste(i, ". Processing", name))
    ica <- CreateSeuratObject(counts = dsets[[i]], project = name)
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
ica.combined <- IntegrateData(anchorset = ica.anchors)
print("Data Integration:")
print(Sys.time() - now)

DefaultAssay(object = ica.combined) <- "integrated"

save(ica.combined, file = "seurat_cca_processed.RData")
print("Saved CCA result to file.")

X <- t(GetAssayData(ica.combined))
writeMM(X, file = "gene_expression.mtx")

write.table(X@Dimnames[[1]], file = "sample_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(X@Dimnames[[2]], file = "feature_names.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
