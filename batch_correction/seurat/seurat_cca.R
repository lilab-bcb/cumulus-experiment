library(Seurat)
library(parallel)
library(R.utils)
library(Matrix)
library(hash)
library(stringi)

seed <- 0
n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores."))

src.obj <- ReadH5AD("../../MantonBM_nonmix_tiny_filter_norm.h5ad")
obj.list <- SplitObject(src.obj, split.by = "Channel")

## Read preset high variable features
df.hvf <- read.csv(file = "../../hvf.txt")
hvf.features <- as.character(unlist(df.hvf['index']))

## Create Seurat objects
obj.list <- lapply(X = obj.list, FUN = function(x) {
	x <- NormalizeData(x)
	x[["RNA"]]@var.features <- hvf.features
	return(x)
})

now <- Sys.time()
ica.anchors <- FindIntegrationAnchors(object.list = obj.list)
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
writeMM(X, file = "matrix.mtx")

write.table(X@Dimnames[[1]], file = "barcodes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(X@Dimnames[[2]], file = "genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
