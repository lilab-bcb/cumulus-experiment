library(Seurat)
library(R.utils)
library(parallel)
library(future)

src_file <- "/data/MantonBM_nonmix_10x.h5"
dst_folder <- "./"
seed <- 0
debug.mode <- TRUE
logfile <- paste0(dst_folder, "seurat_hvg_batch_correction.log")

n.cores <- detectCores()
setOption('mc.cores', n.cores)
write(paste("Use", n.cores, "cores for tSNE."), file = logfile)

options(future.globals.maxSize = 50 * 1024^3)
plan("multiprocess", workers = n.cores)
plan()
write(paste("Use", nbrOfWorkers(), "cores."), file = logfile, append = TRUE)

bm.data <- Read10X_h5(src_file)
n.cells <- dim(bm.data)[2]
bm <- CreateSeuratObject(bm.data, project = "MantonBM", min.cells = floor(n.cells * 0.0005), min.features = 500)
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm <- subset(bm, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mt < 10)

bm.list <- SplitObject(bm)

bm.list <- lapply(X = bm.list, FUN = function(x) {
	x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e5, verbose = debug.mode)
})

now <- Sys.time()
bm.list <- lapply(X = bm.list, FUN = function(x) {
	x <- FindVariableFeatures(x, selection.method = "vst", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000, verbose = debug.mode)
})
logstr.hvg <- Sys.time() - now
write(paste("HVG time:", logstr.hvg, attr(logstr.hvg, "units")), file = logfile, append = TRUE)

plan("multiprocess", workers = 2)
plan()
write(paste("Use", nbrOfWorkers(), "cores."), file = logfile, append = TRUE)

now <- Sys.time()
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, dims = 1:20, verbose = debug.mode)
logstr.anchor <- Sys.time() - now
write(paste("Find Anchors time:", logstr.anchor, attr(logstr.anchor, "units")), file = logfile, append = TRUE)

save(bm.anchors, file = paste0(dst_folder, "seurat_anchor_result.RData"))
rm(bm.list)

now <- Sys.time()
bm.combined <- IntegrateData(anchorset = bm.anchors, dims = 1:20, verbose = debug.mode)
logstr.cca <- Sys.time() - now
write(paste("Integration time:", logstr.cca, attr(logstr.cca, "units")), file = logfile, append = TRUE)

DefaultAssay(bm.combined) <- "integrated"

save(bm.combined, file = paste0(dst_folder, "seurat_corrected_result.RData"))