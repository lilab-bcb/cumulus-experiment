library(Seurat)
library(parallel)
library(future)

seed <- 0
n.cores <- 28
logfile <- "pbmc_seurat_leiden.log"
debug.mode <- FALSE

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan("multiprocess", workers = n.cores)
plan()
write(paste0("Use ", n.cores, " cores."), file = logfile)

load("MantonBM_nonmix.RData")
n.pc <- dim(adata[["pca"]])[2]
adata <- FindNeighbors(adata, k.param = 100, dims = 1:n.pc, nn.method = "annoy", compute.SNN = TRUE, verbose = debug.mode)

print("Finding Clusters using Leiden:")
now <- Sys.time()
adata <- FindClusters(adata, resolution = 1.3, random.seed = seed, algorithm = "leiden", verbose = debug.mode)
logstr.leiden <- Sys.time() - now
write(paste("Leiden time:", logstr.leiden, attr(logstr.leiden, "units")), file = logfile, append = TRUE)

save(adata, file = "seurat_leiden.RData")