sccloud cluster -p 8 --spectral-leiden --umap ../../MantonBM_nonmix_tiny.h5sc ground
sccloud de_analysis -p 8 --labels spectral_leiden_labels --t ground.h5ad ground.de.xlsx
sccloud annotate_cluster ground.h5ad ground.anno.txt
sccloud plot scatter --basis umap --attributes spectral_leiden_labels ground.h5ad ground.umap.pdf