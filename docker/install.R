install.packages('devtools')
install.packages(c('hdf5r', 'leiden', 'uwot', 'ape', 'cowplot', 'fitdistrplus', 'future', 'future.apply', 'ggrepel', 'ggridges', 'ica', 'lmtest', 'metap', 'pbapply', 'plotly', 'png', 'RANN', 'ROCR', 'rsvd', 'Rtsne', 'sctransform', 'SDMTools', 'tsne'))

seurat.url <- "https://cran.r-project.org/src/contrib/Seurat_3.1.0.tar.gz"
install.packages(seurat.url, repos = NULL, type = 'source')