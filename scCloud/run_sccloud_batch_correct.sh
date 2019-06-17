scCloud cluster -p 8 --correct-batch-effect --diffmap-full-speed --run-louvain --run-leiden --run-fitsne --run-umap --run-fle ../MantonBM_nonmix_tiny_10x.h5 MantonBM_nonmix_tiny_sccloud_corrected
scCloud de_analysis -p 8 --labels leiden_labels --mwu --roc MantonBM_nonmix_tiny_sccloud_corrected.h5ad MantonBM_nonmix_tiny_sccloud_corrected.de.xlsx
scCloud annotate_cluster MantonBM_nonmix_tiny_sccloud_corrected.h5ad MantonBM_nonmix_tiny_sccloud_corrected.anno.txt
scCloud plot scatter --basis fitsne --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.fitsne.pdf
scCloud plot scatter --basis umap --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.umap.pdf
scCloud plot scatter --basis fle --attributes leiden_labels,Channel MantonBM_nonmix_tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.fle.pdf