import sccloud as scc
import numpy as np
import pandas as pd
import os, sys
import seaborn as sns
import matplotlib.pyplot as plt
from termcolor import cprint
from sklearn.metrics import adjusted_mutual_info_score

n_cores = os.cpu_count()

src_file = "../MantonBM_nonmix.h5sc"
immune_gene_file = "immune_genes.txt"

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

	cprint("Loading list of immune genes...", "green")
	df_ground = pd.read_csv(immune_gene_file)
	markers = df_ground['gene'].values

	cprint("Processing Seurat HVF selection...", "green")
	adata = scc.read_input(seurat_correct_name + ".h5ad")
	hvf_seurat = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{} genes are marked as HVF in Seurat flavor.".format(hvf_seurat.shape[0]), "yellow")
	marker_seurat = np.intersect1d(hvf_seurat, markers, assume_unique = True)
	cprint("Find {} markers in Seurat HVF selection flavor.".format(marker_seurat.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_seurat}).to_csv("seurat_markers.txt", header = False, index = False)

	cprint("Processing scCloud HVF selection...", "green")
	adata = scc.read_input(sccloud_correct_name + ".h5ad")
	hvf_sccloud = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{} genes are marked as HVF in scCloud flavor.".format(hvf_sccloud.shape[0]), "yellow")
	marker_sccloud = np.intersect1d(hvf_sccloud, markers, assume_unique = True)
	cprint("Find {} markers in scCloud HVF selection flavor.".format(marker_sccloud.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_sccloud}).to_csv("sccloud_markers.txt", header = False, index = False)

def get_mutual_info():

	cprint("Calculating AMI Score on Approximated Leiden Clustering Results between two methods...", "green")

	adata = scc.read_input(seurat_correct_name + ".h5ad")
	bdata = scc.read_input(sccloud_correct_name + ".h5ad")

	mis = adjusted_mutual_info_score(adata.obs['louvain_labels'], bdata.obs['louvain_labels'], average_method = 'arithmetic')
	cprint("AMI = {:.4f}".format(mis))

def plot_figures(file_name):
	if os.system("sccloud plot scatter --basis fitsne --attributes louvain_labels {name}.h5ad {name}.fitsne.pdf".format(name = file_name)):
		sys.exit(1)

	if os.system("sccloud plot scatter --basis fle --attributes louvain_labels {name}.h5ad {name}.fle.pdf".format(name = file_name)):
		sys.exit(1)


if __name__ == '__main__':
	f_list = [f for f in os.listdir('.') if f in [seurat_correct_name + '.h5ad', sccloud_correct_name + '.h5ad']]
	if len(f_list) != 2:
		get_hvf()
		annotate_data(seurat_correct_name)
		annotate_data(sccloud_correct_name)

	compare_markers()
	get_mutual_info()

	plot_figures(seurat_correct_name)
	plot_figures(sccloud_correct_name)

