library(RANN)
library(parallel)
library(R.utils)

n.cores <- detectCores()
setOption('mc.cores', n.cores)
print(paste("Use", n.cores, "cores if possible."))

X <- read.table("x_pca.txt", header = FALSE)
start_seurat <- Sys.time()
my.knn <- nn2(data = X, k = 100, searchtype = "standard", eps = 0)
end_seurat <- Sys.time()
print("Time used for Seurat:")
print(end_seurat - start_seurat)

write.table(my.knn$nn.idx, "seurat_indices.txt", row.names = FALSE, col.names = FALSE)
