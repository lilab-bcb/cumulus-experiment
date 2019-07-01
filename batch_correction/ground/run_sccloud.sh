scCloud cluster -p 8 --run-approximated-leiden --run-fitsne --run-umap ../../MantonBM_nonmix_tiny_10x.h5 ground
scCloud de_analysis -p 8 --labels approx_leiden_labels ground.h5ad ground.de.xlsx
scCloud annotate_cluster ground.h5ad ground.anno.txt
