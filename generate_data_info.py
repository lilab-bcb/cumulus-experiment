import os, sys
import pegasus as pg
import pandas as pd
import numpy as np

from termcolor import cprint

bm_full_data_source = "/data/MantonBM_nonmix.h5sc"
bm_tiny_data_source = "/data/MantonBM_nonmix_tiny.h5sc"
neuron_data_source = "/data/1M_neurons.h5"
pbmc_data_source = '/data/5k_pbmc_v3.h5'

bm_full_out_name = "MantonBM_nonmix_pegasus"

def extract_hvf(input_file, get_pca_knn, batch_correct = True, percent_mito = 10.0):

	fname_prefix = os.path.splitext(input_file)[0]

	# Run correction.
	cprint("Extracting HVGs from {}...".format(input_file), "green")
	adata = pg.read_input(input_file)
	pg.qc_metrics(adata, percent_mito = percent_mito)
	pg.filter_data(adata)
	pg.log_norm(adata)
	pg.highly_variable_features(adata, consider_batch = batch_correct)
	if batch_correct:
		pg.correct_batch(adata, features = "highly_variable_features")

	# Extract highly variable gene set to file.
	cprint("Extracting highly variable features to file...", "green")
	df_hvf = adata.var.loc[adata.var["highly_variable_features"]][["gene_ids"]]
	df_hvf.to_csv(fname_prefix + "_hvf.txt", index_label = "index")

	# PCA and KNN
	if get_pca_knn:
		pg.pca(adata)
		pg.neighbors(adata)
		if batch_correct:
			out_name = fname_prefix + "_corrected"
		else:
			out_name = fname_prefix + "_pca_knn"
		pg.write_output(adata, out_name)

		bdata = adata[:, adata.var['highly_variable_features']].copy()
		if batch_correct:
			out_name_hvf = fname_prefix + "_corrected_hvf"
		else:
			out_name_hvf = fname_prefix + "_pca_knn_hvf"
		pg.write_output(bdata, out_name_hvf)

def preprocess_data(in_file, has_pca, percent_mito = 10.0):

	cprint("Preprocessing {}...".format(in_file), "green")

	# Filtration and Normalization.
	dataset_name = os.path.splitext(in_file)[0].split('/')[-1]
	suffix1 = "_filter_norm"
	output_name1 = os.path.splitext(in_file)[0] + suffix1
	output_hvf_name1 = output_name1 + "_hvf"

	adata = pg.read_input(in_file)
	pg.qc_metrics(adata, percent_mito = percent_mito)
	pg.filter_data(adata)
	pg.log_norm(adata)

	print("Getting Highly Variable Feature Set...")
	df_hvf = pd.read_csv(os.path.splitext(in_file)[0] + "_hvf.txt")
	num_hvf = df_hvf.shape[0]
	df_hvf['highly_variable_features'] = [True] * df_hvf.shape[0]
	df_hvf.set_index('index', inplace = True)
	df_hvf.drop(columns = ['gene_ids'], inplace = True)
	adata.var.drop(columns = ['highly_variable_features'], inplace = True)
	df_join = adata.var.join(df_hvf)
	df_join['highly_variable_features'].fillna(False, inplace = True)
	adata.var['highly_variable_features'] = df_join['highly_variable_features']
	assert np.sum(adata.var['highly_variable_features']) == num_hvf
	adata_hvf = adata[:, adata.var['highly_variable_features']].copy()

	pg.write_output(adata, output_name1)
	pg.write_output(adata_hvf, output_hvf_name1)

	# Furthermore, do PCA.
	suffix2 = "_filter_norm_pca"
	output_name2 = os.path.splitext(in_file)[0] + suffix2
	output_hvf_name2 = output_name2 + "_hvf"

	# Load precalculated PCA
	if has_pca:
		cprint("Set precalculated PCA info...", "green")
		adata.obsm['X_pca'] = np.load("/data/precalculated/{}_uncorrected/pca.npy".format(dataset_name))
		adata.uns['PCs'] = np.load("/data/precalculated/{}_uncorrected/PCs.npy".format(dataset_name))
		adata.uns['pca'] = {}
		adata.uns['pca']['variance'] = np.load("/data/precalculated/{}_uncorrected/pca_variance.npy".format(dataset_name))
		adata.uns['pca']['variance_ratio'] = np.load("/data/precalculated/{}_uncorrected/pca_variance_ratio.npy".format(dataset_name))

		adata_hvf.obsm['X_pca'] = np.load("/data/precalculated/{}_uncorrected/pca.npy".format(dataset_name))
		adata_hvf.uns['PCs'] = np.load("/data/precalculated/{}_uncorrected/PCs.npy".format(dataset_name))
		adata_hvf.uns['pca'] = {}
		adata_hvf.uns['pca']['variance'] = np.load("/data/precalculated/{}_uncorrected/pca_variance.npy".format(dataset_name))
		adata_hvf.uns['pca']['variance_ratio'] = np.load("/data/precalculated/{}_uncorrected/pca_variance_ratio.npy".format(dataset_name))
	else:
		pg.pca(adata)
		pg.pca(adata_hvf)

	pg.write_output(adata, output_name2)
	pg.write_output(adata_hvf, output_hvf_name2)

def annotate_cell_types(label_list):

	adata = pg.read_input(bm_full_out_name + ".h5ad")

	for label, anno_str in label_list:
		anno_dict = {str(i + 1): x for i, x in enumerate(anno_str.split(";"))}
		pg.annotate(adata, 'anno_{}'.format(label), '{}_labels'.format(label), anno_dict)

	pg.write_output(adata, bm_full_out_name)

def run_pegasus_process():
	cprint("Run Pegasus on Bone Marrow dataset with precalculated PCA and Diffusion Map...", "green")

	adata = pg.read_input(bm_full_data_source)
	pg.qc_metrics(adata)
	pg.filter_data(adata)
	pg.log_norm(adata)
	pg.highly_variable_features(adata, consider_batch = True)

	from pegasus.plotting import plot_hvf

	robust_idx = adata.var["robust"].values
	plot_hvf(
		adata.var.loc[robust_idx, "mean"],
 		adata.var.loc[robust_idx, "var"],
		adata.var.loc[robust_idx, "hvf_loess"],
		adata.var.loc[robust_idx, "highly_variable_features"],
		bm_full_out_name + ".hvf.pdf",
	)
	pg.correct_batch(adata, features = "highly_variable_features")

	# Load precalculated PCA
	cprint("Set precalculated PCA info...", "green")
	adata.obsm['X_pca'] = np.load("/data/precalculated/pegasus/pca.npy")
	adata.uns['PCs'] = np.load("/data/precalculated/pegasus/PCs.npy")
	adata.uns['pca'] = {}
	adata.uns['pca']['variance'] = np.load("/data/precalculated/pegasus/pca_variance.npy")
	adata.uns['pca']['variance_ratio'] = np.load("/data/precalculated/pegasus/pca_variance_ratio.npy")

	pg.neighbors(adata)

	# Load precalculated Diffusion Map
	cprint("Set precalculated Diffusion Map info...", "green")
	adata.obsm['X_diffmap'] = np.load("/data/precalculated/pegasus/diffmap.npy")
	adata.uns['diffmap_evals'] = np.load("/data/precalculated/pegasus/diffmap_evals.npy")
	adata.obsm['X_phi'] = np.load("/data/precalculated/pegasus/diffmap_phi.npy")

	pg.write_output(adata, bm_full_out_name)

	del adata

	cprint("Clustering...", "green")
	n_cores = os.cpu_count()
	if os.system("pegasus cluster -p {jobs} --processed --louvain --leiden --spectral-louvain --spectral-leiden --fitsne --tsne --umap --fle --net-tsne --net-umap --net-fle {src}.h5ad null > pegasus.log".format(jobs = n_cores, src = bm_full_out_name)):
		sys.exit(1)

	label_list = ['louvain', 'leiden', 'spectral_louvain', 'spectral_leiden']
	anno_str_louvain = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Naive B cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. Erythrocytes;10. Memory B cells;11. CD14+ Monocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	anno_str_leiden = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. B cells;4. T helper cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. CD14+ Monocytes;10. Erythrocytes;11. Pre B cells;12. HSCs;13. cDCs;14. CD16+ Monocytes;15. Pro B cells;16. pDCs;17. Plasma cells;18. Erythrocytes;19. ;20. Megakaryocytes;21. MSCs"
	anno_str_spectral_louvain = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Cytotoxic T cells;5. Naive B cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. CD14+ Monocytes;10. Memory B cells;11. Erythrocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	anno_str_spectral_leiden = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Naive B cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. Memory B cells;10. CD14+ Monocytes;11. Erythrocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	
	label_list = [('louvain', anno_str_louvain), ('leiden', anno_str_leiden), ('spectral_louvain', anno_str_spectral_louvain), ('spectral_leiden', anno_str_spectral_leiden)]

	cprint("DE Analysis...", "green")
	for label, anno_str in label_list:
		if os.system("pegasus de_analysis -p {jobs} --labels {label}_labels --result-key de_res_{label} --t {src}.h5ad {src}.{label}.de.xlsx > de_{label}.log".format(jobs = n_cores, label = label, src = bm_full_out_name)):
			sys.exit(1)

		if os.system("pegasus annotate_cluster --de-key de_res_{label} {src}.h5ad {src}.{label}.anno.txt > anno_{label}.log".format(label = label, src = bm_full_out_name)):
			sys.exit(1)

	cprint("Annotating...", "green")
	annotate_cell_types(label_list)


if __name__ == '__main__':
	dataset = sys.argv[1]
	assert dataset in ['MantonBM', '1M_neurons', '5k_pbmc']

	if dataset == 'MantonBM':
		extract_hvf(bm_full_data_source, get_pca_knn = True)
		preprocess_data(bm_full_data_source, has_pca = True)

		extract_hvf(bm_tiny_data_source, get_pca_knn = False)
		preprocess_data(bm_tiny_data_source, has_pca = True)

		run_pegasus_process()
	elif dataset == '1M_neurons':
		extract_hvf(neuron_data_source, get_pca_knn = True)
		preprocess_data(neuron_data_source, has_pca = False)
	else:
		extract_hvf(pbmc_data_source, get_pca_knn = True, batch_correct = False, percent_mito = 20)
		preprocess_data(pbmc_data_source, has_pca = False, percent_mito = 20)
