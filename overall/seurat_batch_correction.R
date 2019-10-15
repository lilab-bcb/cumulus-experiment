library(Seurat)
library(future)
library(parallel)

n.cores <- detectCores()
logfile <- "seurat_batch_correction.log"

options(future.globals.maxSize = 50 * 1024^3)  # 50 GB
plan("multiprocess", workers = n.cores)
plan()
write(paste0("Use ", n.cores, " cores."), file = logfile)

## Benchmark Batch Correction. ##
src.obj <- ReadH5AD("/projects/benchmark/MantonBM/MantonBM_nonmix_filter_norm_hvf.h5ad")
VariableFeatures(src.obj) <- GetAssayData(src.obj)@Dimnames[[1]]
obj.list <- SplitObject(src.obj, split.by = "Channel")
now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(obj.list)
logstr.anchors <- Sys.time() - now
write(paste("Finding Anchors time:", logstr.anchors, attr(logstr.anchors, "units")), file = logfile, append = TRUE)

now <- Sys.time()
bm.combined <- IntegrateData(bm.anchors)
logstr.integrate <- Sys.time() - now
write(paste("Integration time:", logstr.integrate, attr(logstr.integrate, "units")), file = logfile, append = TRUE)

save(bm.combined, file = "seurat_cca_result.RData")