scCloud cluster -p 8 --diffmap-full-speed --run-louvain --run-leiden --run-fitsne --run-umap ../MantonBM_nonmix_tiny_10x.h5 MantonBM_nonmix_tiny_sccloud
scCloud de_analysis -p 8 --labels leiden_labels --mwu --roc MantonBM_nonmix_tiny_sccloud.h5ad MantonBM_nonmix_tiny_sccloud.de.xlsx
scCloud annotate_cluster MantonBM_nonmix_tiny_sccloud.h5ad MantonBM_nonmix_tiny_sccloud.anno.txt
scCloud plot scatter --basis fitsne --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud.h5ad tiny_sccloud.fitsne.pdf
scCloud plot scatter --basis umap --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud.h5ad tiny_sccloud.umap.pdf
#scCloud plot scatter --basis fle --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud.h5ad tiny_sccloud.fle.pdf