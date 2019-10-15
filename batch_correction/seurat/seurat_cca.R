library(Seurat)
library(parallel)
library(future)

seed <- 0
n.cores <- detectCores()
logfile <- "seurat_cca.log"

options(future.globals.maxSize = 20 * 1024^3)  # 20 GB
plan("multiprocess", workers = n.cores)
plan()
write(paste("Use", n.cores, "cores."), file = logfile)

src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_tiny_filter_norm.h5ad")
obj.list <- SplitObject(src.obj, split.by = "Channel")
write(paste0("After splitting, there are ", length(obj.list), " batches."), file = logfile, append = TRUE)

## Read preset high variable features
df.hvf <- read.csv(file = "/projects/benchmark/MantonBM/MantonBM_nonmix_hvf.txt")
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
logstr.anchor <- Sys.time() - now
write(paste("Finding Anchors time:", logstr.anchor, attr(logstr.hvgsif, "units")), file = logfile, append = TRUE)

now <- Sys.time()
ica.combined <- IntegrateData(anchorset = ica.anchors)
print("Data Integration:")
logstr.integrate <- Sys.time() - now
write(paste("Integration time:", logstr.integrate, attr(logstr.integrate, "units")), file = logfile, append = TRUE)

DefaultAssay(object = ica.combined) <- "integrated"

save(ica.combined, file = "seurat_cca_processed.RData")
print("Saved CCA result to file.")

X <- t(GetAssayData(ica.combined))
writeMM(X, file = "matrix.mtx")

write.table(X@Dimnames[[1]], file = "barcodes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(X@Dimnames[[2]], file = "genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
