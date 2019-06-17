library(Seurat)

adata <- Read10X_h5("MantonBM_nonmix_10x.h5")

ica <- CreateSeuratObject(counts = adata, project = "ICA", min.cells = 138)
mito.features <- grep(pattern = "^MT-", rownames(ica), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(ica, slot = 'counts')[mito.features,]) / Matrix::colSums(GetAssayData(ica, slot = 'counts'))
ica <- subset(ica, subset = nFeature_RNA >= 500 & nFeature_RNA < 6000 & percent.mito < 0.1)
ica <- NormalizeData(ica, selection.method = "LogNormalize", scale.factor = 1e5)
ica <- FindVariableFeatures(ica, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 7), dispersion.cutoff = c(0.5, Inf))
