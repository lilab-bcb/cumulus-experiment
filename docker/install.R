install.packages('devtools')
install.packages(c('hdf5r', 'leiden', 'uwot'))

seurat.url <- "https://cran.r-project.org/src/contrib/Seurat_3.1.0.tar.gz"
install.packages(seurat.url, repos = NULL, type = 'source')