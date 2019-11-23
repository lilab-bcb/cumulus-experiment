import pegasus as pg
import numpy as np
import pandas as pd
import os, sys
from termcolor import cprint
from sklearn.metrics import adjusted_mutual_info_score

n_cores = os.cpu_count()

src_file = "/data/MantonBM_nonmix.h5sc"
immune_gene_file = "immune_genes.txt"

seurat_correct_name = "MantonBM_nonmix_seurat_hvf_corrected"
pegasus_correct_name = "../MantonBM_nonmix_pegasus"

anno_seurat = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Cytotoxic T cells;5. Naive B cells;6. CD8+ Naive T cells;7. Cytotoxic T cells;8. NK cells;9. CD14+ Monocytes;10. Erythrocytes;11. Memory B cells;12. HSCs;13. Pre B cells;14. cDCs;15. CD16+ Monocytes;16. Erythrocytes;17. Pro B cells;18. pDCs;19. Plasma cells;20. Pre B cells;21. MSCs;22. Cytotoxic T cells"

palettes_seurat = "#c5b0d5,#ff7f0e,#8c564b,#2ca02c,#ff9896,#dbdb8d,#1f77b4,#e377c2,#ffbb78,#9edae5,#aec7e8,#d62728,#98df8a,#9467bd,#c49c94,#ad494a,#f7b6d2,#bcbd22,#17becf,#8c6d31,#000000,#888888"
palettes_pegasus = "#c5b0d5,#ff7f0e,#8c564b,#ff9896,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#9edae5,#aec7e8,#ffbb78,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000"

def get_hvf(processed):

	if processed == 'raw':
		cprint("Computing highly variable genes using Seurat method...", "green")
		if os.system("pegasus cluster -p {jobs} --correct-batch-effect --select-hvf-flavor Seurat --louvain --fitsne {src} {outname}".format(jobs = n_cores, src = src_file, outname = seurat_correct_name)):
			sys.exit(1)
	else:
		cprint("Computing highly variable genes using Seurat method with precalculated PCA...", "green")
		adata =pg.read_input(src_file)
		pg.qc_metrics(adata)
		pg.filter_data(adata)
		pg.log_norm(adata)
		pg.highly_variable_features(adata, consider_batch = True, flavor = 'Seurat')
		pg.correct_batch(adata, features = "highly_variable_features")

		cprint("Set precalculated PCA info...", "green")
		adata.obsm['X_pca'] = np.load("/data/precalculated/seurat_hvf/pca.npy")
		adata.uns['PCs'] = np.load("/data/precalculated/seurat_hvf/PCs.npy")
		adata.uns['pca'] = {}
		adata.uns['pca']['variance'] = np.load("/data/precalculated/seurat_hvf/pca_variance.npy")
		adata.uns['pca']['variance_ratio'] = np.load("/data/precalculated/seurat_hvf/pca_variance_ratio.npy")

		pg.neighbors(adata)

		pg.write_output(adata, seurat_correct_name)
		del adata

		cprint("Clustering...", "green")
		if os.system("pegasus cluster -p {jobs} --processed --louvain --fitsne {src} {outname}".format(jobs = n_cores, src = seurat_correct_name + '.h5ad', outname = seurat_correct_name)):
			sys.exit(1)


def annotate_data(processed):
	cprint("Annotating Cells for {name}...".format(name = seurat_correct_name), "green")

	if os.system("pegasus de_analysis -p {jobs} --labels louvain_labels --t {name}.h5ad {name}.de.xlsx".format(jobs = n_cores, name = seurat_correct_name)):
		sys.exit(1)

	if os.system("pegasus annotate_cluster {name}.h5ad {name}.anno.txt".format(name = seurat_correct_name)):
		sys.exit(1)

	if processed == 'processed':
		adata = pg.read_input(seurat_correct_name + '.h5ad')
		anno_dict = {str(i + 1): x for i, x in enumerate(anno_seurat.split(";"))}
		pg.annotate(adata, 'anno', 'louvain_labels', anno_dict)
		pg.write_output(adata, seurat_correct_name)


def compare_markers():

	cprint("Loading list of immune genes...", "green")
	df_ground = pd.read_csv(immune_gene_file)
	markers = df_ground['gene'].values

	cprint("Processing Seurat HVF selection...", "green")
	adata = pg.read_input(seurat_correct_name + ".h5ad")
	hvf_seurat = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{} genes are marked as HVF in Seurat flavor.".format(hvf_seurat.shape[0]), "yellow")
	marker_seurat = np.intersect1d(hvf_seurat, markers, assume_unique = True)
	cprint("Find {} markers in Seurat HVF selection flavor.".format(marker_seurat.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_seurat}).to_csv("seurat_markers.txt", header = False, index = False)

	cprint("Processing Pegasus HVF selection...", "green")
	adata = pg.read_input(pegasus_correct_name + ".h5ad")
	hvf_pegasus = adata.var.loc[adata.var['highly_variable_features']].index.values
	cprint("{} genes are marked as HVF in Pegasus flavor.".format(hvf_pegasus.shape[0]), "yellow")
	marker_pegasus = np.intersect1d(hvf_pegasus, markers, assume_unique = True)
	cprint("Find {} markers in Pegasus HVF selection flavor.".format(marker_pegasus.shape[0]), "yellow")
	pd.DataFrame({'gene': marker_pegasus}).to_csv("pegasus_markers.txt", header = False, index = False)

	cprint("Processing common markers that both flavors find...", "green")
	common_markers = np.intersect1d(marker_seurat, marker_pegasus, assume_unique = True)
	cprint("{} markers are commonly found by both flavors.".format(common_markers.shape[0]), "yellow")
	pd.DataFrame({'gene': common_markers}).to_csv("common_markers.txt", header = False, index = False)

	cprint("Processing Seurat flavor specific markers...", "green")
	seurat_specific = np.setdiff1d(marker_seurat, common_markers, assume_unique = True)
	cprint("{} markers are Seurat flavor specific.".format(seurat_specific.shape[0]), "yellow")
	pd.DataFrame({'gene': seurat_specific}).to_csv("seurat_specific.txt", header = False, index = False)

	cprint("Processing Pegasus flavor specific markers...", "green")
	pegasus_specific = np.setdiff1d(marker_pegasus, common_markers, assume_unique = True)
	cprint("{} markers are Pegasus flavor specific.".format(pegasus_specific.shape[0]), "yellow")
	pd.DataFrame({'gene': pegasus_specific}).to_csv("pegasus_specific.txt", header = False, index = False)

def get_mutual_info():

	cprint("Calculating AMI Score on Approximated Leiden Clustering Results between two methods...", "green")

	adata = pg.read_input(seurat_correct_name + ".h5ad")
	bdata = pg.read_input(pegasus_correct_name + ".h5ad")

	ami = adjusted_mutual_info_score(adata.obs['louvain_labels'], bdata.obs['louvain_labels'], average_method = 'arithmetic')
	cprint("AMI = {:.4f}".format(ami), "yellow")

def plot_figures(processed):

	if processed == 'processed':
		if os.system('pegasus plot scatter --basis fitsne --attributes anno --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/Figure_S1C_left.pdf'.format(src = seurat_correct_name, palettes = palettes_seurat)):
			sys.exit(1)
	else:
		if os.system('pegasus plot scatter --basis fitsne --attributes louvain_labels {src}.h5ad /output/Figure_S1C_left.pdf'.format(src = seurat_correct_name)):
			sys.exit(1)

	if os.system('pegasus plot scatter --basis fitsne --attributes anno_louvain --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/Figure_S1C_right.pdf'.format(src = pegasus_correct_name, palettes = palettes_pegasus)):
		sys.exit(1)

	if os.system("cp {src}.hvf.pdf /output/Figure_S1A.pdf".format(src = pegasus_correct_name)):
		sys.exit(1)


if __name__ == '__main__':
	processed = sys.argv[1]
	assert processed in ['processed', 'raw']
	f_list = [f for f in os.listdir('.') if f in [seurat_correct_name + '.h5ad']]
	if len(f_list) != 1:
		get_hvf(processed)
		annotate_data(processed)

	compare_markers()
	get_mutual_info()

	plot_figures(processed)

