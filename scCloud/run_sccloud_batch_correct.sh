scCloud cluster -p 8 --correct-batch-effect --diffmap-full-speed --run-louvain --run-leiden --run-fitsne --run-umap ../MantonBM_nonmix_tiny_10x_filtered.h5ad MantonBM_nonmix_tiny_filtered_corrected
#scCloud de_analysis -p 8 --labels leiden_labels --mwu --roc MantonBM_nonmix_tiny_filtered_corrected.h5ad MantonBM_nonmix_tiny_filtered_corrected.de.xlsx
#scCloud annotate_cluster MantonBM_nonmix_tiny_filtered_corrected.h5ad MantonBM_nonmix_tiny_filtered_corrected.anno.txt
#scCloud plot scatter --basis fitsne --attributes leiden_labels,Channel MantonBM_nonmix_tiny_filtered_corrected.h5ad tiny_filtered_corrected.fitsne.pdf
#scCloud plot scatter --basis umap --attributes leiden_labels,Channel MantonBM_nonmix_tiny_filtered_corrected.h5ad tiny_filtered_corrected.umap.pdf