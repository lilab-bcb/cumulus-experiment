import sccloud as scc
import numpy as np
import pandas as pd
import os, sys
import seaborn as sns
import matplotlib.pyplot as plt
from termcolor import cprint
from sklearn.metrics import adjusted_mutual_info_score

n_cores = os.cpu_count()

markers = np.array(['CD3D', 'CD3E', 'CD3G', 'TRAC', 'CD4', 'CD8A', 'CD8B', 'RTKN2', 'TIGIT', 'FOXP3', 'IL2RA', 'CCR7', 'SELL', 'IL7R', 'TCF7', 'CD27', 'CD19', 'MS4A1', 'CD79A', 'CD79B', 'MZB1', 'HLA-DRA', 'HLA-DRB1', 'CD34', 'IGLL1', 'NCAM1', 'NKG7', 'KLRB1', 'KLRD1', 'KLRF1', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRC4', 'FCGR3A', 'ITGAL', 'ITGAM', 'VCAN', 'FCN1', 'S100A8', 'S100A9', 'CD14', 'ASAH1', 'MS4A7', 'IFITM2', 'IFITM3', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'FCER1A', 'CLEC10A', 'CD1C', 'THBD', 'JCHAIN', 'LILRA4', 'GZMB', 'IL3RA', 'SERPINF1', 'ITM2C', 'IRF7', 'CD38', 'XBP1', 'SLAMF7', 'TNFRSF17', 'TNFRSF13B', 'IGHA1', 'IGHG1', 'KIT', 'CD59', 'THY1', 'SOX4', 'GYPA', 'HBB', 'HBA1', 'TFRC', 'ITGA4', 'ANK1', 'ICAM4', 'BCAM', 'SLC4A1', 'ACKR1', 'PF4', 'PPBP', 'GP5', 'CXCR4', 'SLAMF1', 'MPL', 'ITGA2B', 'FUT4', 'MPO', 'CSF3R', 'FCGR3B', 'CEACAM8'])

src_file = "../MantonBM_nonmix.h5sc"

seurat_correct_name = "MantonBM_nonmix_seurat_corrected"
sccloud_correct_name = "MantonBM_nonmix_sccloud_corrected"


def get_hvf():
	cprint("Computing highly variable genes using Seurat method...", "green")
	if os.system("sccloud cluster -p {jobs} --correct-batch-effect --select-hvf-flavor Seurat --louvain --fitsne --fle {src} {outname}".format(jobs = n_cores, src = src_file, outname = seurat_correct_name)):
		sys.exit(1)

	cprint("Computing highly variable genes using scCloud new method...", "green")
	if os.system("sccloud cluster -p {jobs} --plot-hvf --correct-batch-effect --louvain --fitsne --fle {src} {outname}".format(jobs = n_cores, src = src_file, outname = sccloud_correct_name)):
		sys.exit(1)


def annotate_data(file_name):
	cprint("Annotating Cells for {name}...".format(name = file_name), "green")

	if os.system("sccloud de_analysis -p {jobs} --labels louvain_labels --t {name}.h5ad {name}.de.xlsx".format(jobs = n_cores, name = file_name)):
		sys.exit(1)

	if os.system("sccloud annotate_cluster {name}.h5ad {name}.anno.txt".format(name = file_name)):
		sys.exit(1)


def compare_markers():

	cprint("Processing Seurat HVF selection...", "green")

	adata = scc.read_input(seurat_correct_name + ".h5ad", mode = 'a')
	hvf_seurat = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{num} genes are marked as HVF in Seurat method.".format(num = hvf_seurat.shape[0]), "yellow")
	pd.DataFrame({'gene': hvf_seurat}).to_csv("seurat_hvf.txt", header = False, index = False)
	marker_seurat = np.intersect1d(hvf_seurat, markers, assume_unique = True)
	cprint("Find {num} markers in Seurat HVF selection method.".format(num = marker_seurat.shape[0]), "yellow")

	cprint("Processing scCloud HVF selection...", "green")
	adata = scc.read_input(sccloud_correct_name + ".h5ad", mode = 'a')
	hvf_sccloud = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{num} genes are marked as HVF in scCloud method.".format(num = hvf_sccloud.shape[0]), "yellow")
	pd.DataFrame({'gene': hvf_sccloud}).to_csv("sccloud_hvf.txt", header = False, index = False)
	marker_sccloud = np.intersect1d(hvf_sccloud, markers, assume_unique = True)
	cprint("Find {num} markers in scCloud HVF selection method.".format(num = marker_sccloud.shape[0]), "yellow")

	common_markers = np.intersect1d(marker_seurat, marker_sccloud, assume_unique = True)
	cprint("Among them, {num} markers are found by both methods.".format(num = common_markers.shape[0]), "yellow")
	pd.DataFrame({'gene': common_markers}).to_csv("common_markers.txt", header = False, index = False)

	# Compute set difference.
	hvf_seurat_only = np.setdiff1d(marker_seurat, common_markers)
	cprint("{num} markers are found by Seurat method only.".format(num = hvf_seurat_only.shape[0]), "yellow")
	pd.DataFrame({'gene': hvf_seurat_only}).to_csv("seurat_only_marker.txt", header = False, index = False)

	hvf_sccloud_only = np.setdiff1d(marker_sccloud, common_markers)
	cprint("{num} markers are found by scCloud method only.".format(num = hvf_sccloud_only.shape[0]), "yellow")
	pd.DataFrame({'gene': hvf_sccloud_only}).to_csv("sccloud_only_marker.txt", header = False, index = False)

	# Compute markers not covered by either method.
	marker_uncovered = np.setdiff1d(markers, np.union1d(marker_seurat, marker_sccloud))
	cprint("{} markers are not covered by either method.".format(marker_uncovered.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_uncovered}).to_csv("marker_uncovered.txt", header = False, index = False)

	# Write all markers to csv.
	cprint("{} markers in total.".format(markers.shape[0]), "yellow")
	pd.DataFrame({'gene': markers}).to_csv("all_markers.txt", header = False, index = False)

def get_mutual_info():

	cprint("Calculating AMI Score on Approximated Leiden Clustering Results between two methods...", "green")

	adata = scc.read_input(seurat_correct_name + ".h5ad", mode = 'a')
	bdata = scc.read_input(sccloud_correct_name + ".h5ad", mode = 'a')

	mis = adjusted_mutual_info_score(adata.obs['louvain_labels'], bdata.obs['louvain_labels'], average_method = 'arithmetic')
	cprint("AMI = {:.4f}".format(mis))

def plot_figures():
	if os.system("sccloud plot scatter --basis fle --attributes louvain_labels {name}.h5ad {name}.fle.pdf".format(name = seurat_correct_name)):
		sys.exit(1)

	if os.system("sccloud plot scatter --basis fle --attributes louvain_labels {name}.h5ad {name}.fle.pdf".format(name = sccloud_correct_name)):
		sys.exit(1)


def plot_diffusion_coeffs():

	# For alpha = 0.5
	adata = scc.read_input(sccloud_correct_name + '.h5ad', mode = 'a')
	S = adata.uns['diffmap_evals']
	alpha = 0.5
	weights = S / (1 - alpha * S)

	fig, ax = plt.subplots()
	fig.set_size_inches(18, 8)
	sns.barplot(x = list(range(1, weights.shape[0] + 1)), y = weights, color = 'salmon', ax = ax)
	ax.set_xlabel("Diffusion Components", fontsize = 15)
	ax.set_ylabel("Diffusion Coefficients", fontsize = 15)
	ax.set(ylim = (0, 2))
	fig.savefig("new_corrected_randomize_alpha_0.5.dist.pdf")
	plt.close()

	# For alpha = 1.0
	bdata = scc.read_input(sccloud_correct_alpha_one_name + '.h5ad', mode = 'a')
	S = bdata.uns['diffmap_evals']
	alpha = 1.0
	weights = S / (1 - alpha * S)

	fig, ax = plt.subplots()
	fig.set_size_inches(18, 8)
	sns.barplot(x = list(range(1, weights.shape[0] + 1)), y = weights, color = 'royalblue', ax = ax)
	ax.set_xlabel("Diffusion Components", fontsize = 15)
	ax.set_ylabel("Diffusion Coefficients", fontsize = 15)
	ax.set(ylim = (0, 350))
	fig.savefig("new_corrected_randomize_alpha_1.0.dist.pdf")
	plt.close()


if __name__ == '__main__':
	f_list = [f for f in os.listdir('.') if f in [seurat_correct_name + '.h5ad', sccloud_correct_name + '.h5ad']]
	if len(f_list) != 2:
		get_hvf()
		annotate_data(seurat_correct_name)
		annotate_data(sccloud_correct_name)

	compare_markers()
	get_mutual_info()

	plot_figures()

	#plot_diffusion_coeffs()
