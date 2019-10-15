library(Seurat)
library(future)
library(parallel)

n.cores <- detectCores()
logfile <- "seurat_batch_correction.log"

options(future.globals.maxSize = 50 * 1024^3)  # 50 GB
plan("multiprocess", workers = n.cores)
plan()

## Benchmark Batch Correction. ##
src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_filter_norm_hvf.h5ad")
VariableFeatures(src.obj) <- GetAssayData(src.obj)@Dimnames[[1]]
obj.list <- SplitObject(src.obj, split.by = "Channel")
now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(obj.list)
logstr.knn <- Sys.time() - now
write(paste("KNN time:", logstr.knn, attr(logstr.knn, "units")), file = logfile, append = TRUE)

bm.combined <- IntegrateData(bm.anchors)
print("Batch Correction:")
print(Sys.time() - now)