sccloud cluster -p 8 --louvain --umap ../../MantonBM_nonmix_tiny.h5sc ground
sccloud de_analysis -p 8 --labels louvain_labels --t ground.h5ad ground.de.xlsx
sccloud annotate_cluster ground.h5ad ground.anno.txt
sccloud plot scatter --basis umap --attributes louvain_labels,Channel ground.h5ad ground.umap.pdf