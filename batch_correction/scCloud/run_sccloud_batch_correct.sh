scCloud cluster -p 8 --correct-batch-effect --diffmap-full-speed --run-louvain --run-leiden --run-approximated-louvain --run-approximated-leiden --run-fitsne --run-umap ../../MantonBM_nonmix_tiny_10x.h5 tiny_sccloud_corrected
#scCloud de_analysis -p 8 --labels leiden_labels --mwu --roc tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.de.xlsx
#scCloud find_markers -p 8 --labels leiden_labels --remove-ribo tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.markers.xlsx
#scCloud annotate_cluster tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.anno.txt
#scCloud plot scatter --basis fitsne --attributes leiden_labels,Channel tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.fitsne.pdf
#scCloud plot scatter --basis umap --attributes leiden_labels,Channel tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.umap.pdf