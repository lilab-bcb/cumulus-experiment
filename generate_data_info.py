import os, sys
import pegasus as pg
import pandas as pd
from termcolor import cprint

bm_full_data_source = "/data/MantonBM_nonmix.h5sc"
bm_tiny_data_source = "/data/MantonBM_nonmix_tiny.h5sc"
neuron_data_source = "/data/1M_neurons.h5"

bm_full_out_name = "MantonBM_nonmix_pegasus"

def extract_hvf(input_file):

	fname_prefix = os.path.splitext(input_file)[0]

	# Run correction.
	cprint("Extracting HVGs from {}...".format(input_file), "green")
	adata = pg.read_input(input_file)
	pg.qc_metrics(adata)
	pg.filter_data(adata)
	pg.log_norm(adata)
	pg.highly_variable_features(adata, consider_batch = True)
	pg.correct_batch(adata, features = "highly_variable_features")

	# Extract highly variable gene set to file.
	cprint("Extracting highly variable features to file...", "green")
	df_hvf = adata.var.loc[adata.var["highly_variable_features"]][["gene_ids"]]
	df_hvf.to_csv(fname_prefix + "_hvf.txt", index_label = "index")

	# PCA and KNN
	pg.pca(adata)
	pg.neighbors(adata)
	pg.write_output(adata, fname_prefix + "_corrected")

	bdata = adata[:, adata.var['highly_variable_features']].copy()
	pg.write_output(bdata, fname_prefix + "_corrected_hvf")

def preprocess_data(in_file):

	cprint("Preprocessing {}...".format(in_file), "green")

	# Filtration and Normalization.
	suffix1 = "_filter_norm"
	output_name1 = os.path.splitext(in_file)[0] + suffix1
	output_hvf_name1 = output_name1 + "_hvf"

	adata = pg.read_input(in_file)
	pg.qc_metrics(adata)
	pg.filter_data(adata)
	pg.log_norm(adata)

	print("Getting Highly Variable Feature Set...")
	in_name = os.path.splitext(in_file)[0]
	if in_name.split('_')[-1] == "tiny":
		df_hvf = pd.read_csv(in_name[0:-5] + "_hvf.txt")
	else:
		df_hvf = pd.read_csv(os.path.splitext(in_file)[0] + "_hvf.txt")
	df_hvf['highly_variable_features'] = [True] * df_hvf.shape[0]
	df_hvf.set_index('index', inplace = True)
	df_hvf.drop(columns = ['gene_ids'], inplace = True)
	adata.var.drop(columns = ['highly_variable_features'], inplace = True)
	df_join = adata.var.join(df_hvf)
	df_join['highly_variable_features'].fillna(False, inplace = True)
	adata.var['highly_variable_features'] = df_join['highly_variable_features']
	adata_hvf = adata[:, adata.var['highly_variable_features']].copy()

	pg.write_output(adata, output_name1)
	pg.write_output(adata_hvf, output_hvf_name1)

	# Furthermore, do PCA.
	suffix2 = "_filter_norm_pca"
	output_name2 = os.path.splitext(in_file)[0] + suffix2
	output_hvf_name2 = output_name2 + "_hvf"

	pg.pca(adata)
	pg.write_output(adata, output_name2)

	pg.pca(adata_hvf)
	pg.write_output(adata_hvf, output_hvf_name2)

def run_pegasus_process():
	cprint("Run pegasus on Bone Marrow dataset...", "green")
	n_cores = os.cpu_count()
	if os.system("pegasus cluster -p {jobs} --plot-hvf --correct-batch-effect --diffmap --louvain --leiden --spectral-louvain --spectral-leiden --fitsne --tsne --umap --fle --net-tsne --net-umap --net-fle {src} {dst} > pegasus.log".format(jobs = n_cores, src = bm_full_data_source, dst = bm_full_out_name)):
		sys.exit(1)

	label_list = ['louvain', 'leiden', 'spectral_louvain', 'spectral_leiden']
	anno_str_louvain = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Naive B cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. Erythrocytes;10. Memory B cells;11. CD14+ Monocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	anno_str_leiden = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. B cells;4. T helper cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. CD14+ Monocytes;10. Erythrocytes;11. Pre B cells;12. HSCs;13. cDCs;14. CD16+ Monocytes;15. Pro B cells;16. pDCs;17. Plasma cells;18. Erythrocytes;19. ;20. Megakaryocytes;21. MSCs"
	anno_str_spectral_louvain = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Cytotoxic T cells;5. Naive B cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. CD14+ Monocytes;10. Memory B cells;11. Erythrocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	anno_str_spectral_leiden = "1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Naive B cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. Memory B cells;10. CD14+ Monocytes;11. Erythrocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"
	
	label_list = [('louvain', anno_str_louvain), ('leiden', anno_str_leiden), ('spectral_louvain', anno_str_spectral_louvain), ('spectral_leiden', anno_str_spectral_leiden)]

	for label, anno_str in label_list:
		if os.system("pegasus de_analysis -p {jobs} --labels {label}_labels --result-key de_res_{label} --t {src}.h5ad {src}.{label}.de.xlsx > de_{label}.log".format(jobs = n_cores, label = label, src = bm_full_out_name)):
			sys.exit(1)

		if os.system("pegasus annotate_cluster --de-key de_res_{label} {src}.h5ad {src}.{label}.anno.txt > anno_{label}.log".format(label = label, src = bm_full_out_name)):
			sys.exit(1)

		if os.system('pegasus annotate_cluster --annotation "anno_{label}:{label}_labels:{anno}" {src}.h5ad'.format(label = label, anno = anno_str, src = bm_full_out_name)):
			sys.exit(1)

if __name__ == '__main__':
	dataset = sys.argv[1]
	assert dataset in ['MantonBM', '1M_neurons']

	if dataset == 'MantonBM':
		extract_hvf(bm_full_data_source)
		preprocess_data(bm_full_data_source)
		preprocess_data(bm_tiny_data_source)
		run_pegasus_process()
	else:
		extract_hvf(neuron_data_source)
		preprocess_data(neuron_data_source)
