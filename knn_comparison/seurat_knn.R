library(RANN)
library(parallel)
library(future.apply)

source("/opt/software/seurat-3.1.0/R/clustering.R")

n.cores <- detectCores()
plan("multiprocess", workers = n.cores)
plan()

X <- read.table("x_pca.txt", header = FALSE)

## For nn2 ##
#start_ann <- Sys.time()
#my.knn1 <- nn2(data = X, k = 100, searchtype = "standard", eps = 0)
#end_ann <- Sys.time()
#print("Time used for Seurat ANN:")
#print(end_ann - start_ann)
#write.table(my.knn1$nn.idx, "seurat_indices_ann.txt", row.names = FALSE, col.names = FALSE)

## For AnnoyNN ##
Y <- as.matrix(X)
start_annoy <- Sys.time()
my.knn2 <- AnnoyNN(data = Y, k = 100)
logstr.annoy <- Sys.time() - start_annoy
print(paste("Time used for Seurat AnnoyNN:", logstr.annoy, attr(logstr.annoy, "units")))
write.table(my.knn2$nn.idx, "seurat_indices_annoy.txt", row.names = FALSE, col.names = FALSE)

save(my.knn2, file = "seurat_knn_result.RData")