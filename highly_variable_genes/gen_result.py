import scCloud
import numpy as np
import pandas as pd
import os, sys
from termcolor import cprint
from sklearn.metrics import adjusted_mutual_info_score

markers = np.array(['CD3D', 'CD3E', 'CD3G', 'TRAC', 'CD4', 'CD8A', 'CD8B', 'RTKN2', 'TIGIT', 'FOXP3', 'IL2RA', 'CCR7', 'SELL', 'IL7R', 'TCF7', 'CD27', 'CD19', 'MS4A1', 'CD79A', 'CD79B', 'MZB1', 'HLA-DRA', 'HLA-DRB1', 'CD34', 'IGLL1', 'NCAM1', 'NKG7', 'KLRB1', 'KLRD1', 'KLRF1', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRC4', 'FCGR3A', 'ITGAL', 'ITGAM', 'VCAN', 'FCN1', 'S100A8', 'S100A9', 'CD14', 'ASAH1', 'MS4A7', 'IFITM2', 'IFITM3', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'FCER1A', 'CLEC10A', 'CD1C', 'THBD', 'JCHAIN', 'LILRA4', 'GZMB', 'IL3RA', 'SERPINF1', 'ITM2C', 'IRF7', 'CD38', 'XBP1', 'SLAMF7', 'TNFRSF17', 'TNFRSF13B', 'IGHA1', 'IGHG1', 'KIT', 'CD59', 'THY1', 'SOX4', 'GYPA', 'HBB', 'HBA1', 'TFRC', 'ITGA4', 'ANK1', 'ICAM4', 'BCAM', 'SLC4A1', 'ACKR1', 'PF4', 'PPBP', 'GP5', 'CXCR4', 'SLAMF1', 'MPL', 'ITGA2B', 'FUT4', 'MPO', 'CSF3R', 'FCGR3B', 'CEACAM8'])

seurat_correct_name = "MantonBM_nonmix_seurat_corrected"
sccloud_correct_name = "MantonBM_nonmix_new_corrected"


def get_hvg():
	cprint("Computing highly variable genes using Seurat method...", "green")
	if os.system("scCloud cluster -p 8 --correct-batch-effect --select-hvg-flavor Seurat --run-leiden --run-approximated-leiden --run-louvain --run-approximated-louvain --run-tsne --run-fitsne --run-umap --run-fle ../MantonBM_nonmix_10x.h5 {outname}".format(outname = seurat_correct_name)):
		sys.exit(1)

	cprint("Computing highly variable genes using scCloud new method...", "green")
	if os.system("scCloud cluster -p 8 --plot-hvg --correct-batch-effect --run-leiden --run-approximated-leiden --run-louvain --run-approximated-louvain --run-tsne --run-fitsne --run-umap --run-fle ../MantonBM_nonmix_10x.h5 {outname}".format(outname = sccloud_correct_name)):
		sys.exit(1)


def annotate_data(file_name):
	cprint("Annotating Cells for {name}...".format(name = file_name), "green")

	if os.system("scCloud de_analysis -p 8 --labels leiden_labels --fisher --mwu --roc {name}.h5ad {name}.de.xlsx".format(name = file_name)):
		sys.exit(1)

	if os.system("scCloud annotate_cluster {name}.h5ad {name}.anno.txt".format(name = file_name)):
		sys.exit(1)


def compare_markers():

	cprint("Processing Seurat HVG selection...", "green")

	adata = scCloud.tools.read_input(seurat_correct_name + ".h5ad", mode = 'a')
	hvg_seurat = adata.var.loc[adata.var['highly_variable_genes']].index.values
	cprint("{num} genes are marked as HVG in Seurat method.".format(num = hvg_seurat.shape[0]), "yellow")
	pd.DataFrame({'gene': hvg_seurat}).to_csv("seurat_hvg.txt", header = False, index = False)
	marker_seurat = np.intersect1d(hvg_seurat, markers, assume_unique = True)
	cprint("Find {num} markers in Seurat HVG selection method.".format(num = marker_seurat.shape[0]), "yellow")

	cprint("Processing scCloud HVG selection...", "green")
	adata = scCloud.tools.read_input(sccloud_correct_name + ".h5ad", mode = 'a')
	hvg_sccloud = adata.var.loc[adata.var['highly_variable_genes']].index.values
	cprint("{num} genes are marked as HVG in scCloud method.".format(num = hvg_sccloud.shape[0]), "yellow")
	pd.DataFrame({'gene': hvg_sccloud}).to_csv("sccloud_hvg.txt", header = False, index = False)
	marker_sccloud = np.intersect1d(hvg_sccloud, markers, assume_unique = True)
	cprint("Find {num} markers in scCloud HVG selection method.".format(num = marker_sccloud.shape[0]), "yellow")

	common_markers = np.intersect1d(marker_seurat, marker_sccloud, assume_unique = True)
	cprint("Among them, {num} markers are found by both methods.".format(num = common_markers.shape[0]), "yellow")
	pd.DataFrame({'gene': common_markers}).to_csv("common_markers.txt", header = False, index = False)

	# Compute set difference.
	hvg_seurat_only = np.setdiff1d(marker_seurat, common_markers)
	cprint("{num} markers are found by Seurat method only.".format(num = hvg_seurat_only.shape[0]), "yellow")
	pd.DataFrame({'gene': hvg_seurat_only}).to_csv("seurat_only_marker.txt", header = False, index = False)

	hvg_sccloud_only = np.setdiff1d(marker_sccloud, common_markers)
	cprint("{num} markers are found by scCloud method only.".format(num = hvg_sccloud_only.shape[0]), "yellow")
	pd.DataFrame({'gene': hvg_sccloud_only}).to_csv("sccloud_only_marker.txt", header = False, index = False)

	# Compute markers not covered by either method.
	marker_uncovered = np.setdiff1d(markers, np.union1d(marker_seurat, marker_sccloud))
	cprint("{} markers are not covered by either method.".format(marker_uncovered.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_uncovered}).to_csv("marker_uncovered.txt", header = False, index = False)

	# Write all markers to csv.
	cprint("{} markers in total.".format(markers.shape[0]), "yellow")
	pd.DataFrame({'gene': markers}).to_csv("all_markers.txt", header = False, index = False)

def get_mutual_info():

	cprint("Calculating AMI Score on Approximated Leiden Clustering Results between two methods...", "green")

	adata = scCloud.tools.read_input(seurat_correct_name + ".h5ad", mode = 'a')
	bdata = scCloud.tools.read_input(sccloud_correct_name + ".h5ad", mode = 'a')

	mis = adjusted_mutual_info_score(adata.obs['approx_leiden_labels'], bdata.obs['approx_leiden_labels'], average_method = 'arithmetic')
	cprint("AMI = {:.4f}".format(mis))

def plot_figures(file_name):
	if os.system("scCloud plot scatter --basis tsne --attributes approx_leiden_labels,Individual {name}.h5ad {name}.tsne.pdf".format(name = file_name)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis fitsne --attributes approx_leiden_labels,Individual {name}.h5ad {name}.fitsne.pdf".format(name = file_name)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis umap --attributes approx_leiden_labels,Individual {name}.h5ad {name}.umap.pdf".format(name = file_name)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis fle --attributes approx_leiden_labels,Individual {name}.h5ad {name}.fle.pdf".format(name = file_name)):
		sys.exit(1)


if __name__ == '__main__':
	f_list = [f for f in os.listdir('.') if f in ["{}.h5ad".format(seurat_correct_name), "{}.h5ad".format(sccloud_correct_name)]]
	if len(f_list) != 2:
		get_hvg()
		annotate_data(seurat_correct_name)
		annotate_data(sccloud_correct_name)

	compare_markers()
	get_mutual_info()

	#plot_figures(seurat_correct_name)
	#plot_figures(sccloud_correct_name)
