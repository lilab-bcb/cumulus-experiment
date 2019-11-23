import os, sys
import pegasus as pg
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import time
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"
data_dst = "spectral_result"
spectral_label = "spectral_labels"

anno_str = "1. CD14+ Monocytes;2. CD4+ Naive T cells;3. B cells;4. Erythrocytes;5. Plasma cells;6. Megakaryocytes;7. NK cells;8. Pre B cells;9. pDCs;10. Erythrocytes;11. Pro B cells;12. CD16+ Monocytes;13. MSCs;14. CD14+ Monocytes;15. HSCs;16. Pre B cells;17. HSCs;18. Erythrocytes;19. cDCs;20. Cytotoxic T cells;21. Erythrocytes"

palettes_spectral = "#ffbb78,#c5b0d5,#ff9896,#9edae5,#17becf,#8c6d31,#e377c2,#98df8a,#bcbd22,#125b69,#f7b6d2,#c49c94,#000000,#ff7f0e,#d62728,#1cda9b,#f705ec,#ad494a,#9467bd,#1f77b4,#250bca"

def run_spectral(data, rep_key, n_clusters, K = 100, n_jobs = 1, random_state = 0):

	X = data.obsm[rep_key].astype('float64')
	km = KMeans(n_clusters = n_clusters, n_jobs = n_jobs, random_state = random_state)
	km.fit(X)

	data.obs[spectral_label] = pd.Categorical((km.labels_.astype('int') + 1).astype('str'))


if __name__ == '__main__':

	if (data_dst + '.h5ad') not in os.listdir('.'):
		adata = pg.read_input(data_src + '.h5ad')
		start_spec = time.time()
		run_spectral(adata, 'X_diffmap', n_clusters = 21, n_jobs = n_cores)
		end_spec = time.time()
		cprint("Time for spectral clustering = {} s.".format(end_spec - start_spec), "yellow")
		pg.fitsne(adata, rep = 'pca', n_jobs = n_cores)
		pg.write_output(adata, data_dst)

	if os.system("pegasus de_analysis -p {jobs} --labels {label} --t {name}.h5ad {name}.de.xlsx".format(jobs = n_cores, label = spectral_label, name = data_dst)):
		sys.exit(1)

	if os.system("pegasus annotate_cluster {name}.h5ad {name}.anno.txt".format(name = data_dst)):
		sys.exit(1)

	adata = pg.read_input(data_dst + '.h5ad')
	anno_dict = {str(i + 1): x for i, x in enumerate(anno_str.split(";"))}
	pg.annotate(adata, 'anno', spectral_label, anno_dict)
	pg.write_output(adata, data_dst)

	if os.system('pegasus plot scatter --basis fitsne --attributes anno --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/Figure_S5A_left.pdf'.format(palettes = palettes_spectral, src = data_dst)):
		sys.exit(1)

	adata = pg.read_input(data_dst + '.h5ad')
	cprint("Calculating AMI...", "green")
	ami_spectral_vs_louvain = adjusted_mutual_info_score(adata.obs[spectral_label], adata.obs['louvain_labels'], average_method = 'arithmetic')
	cprint("AMI between {label} and louvain_labels: {ami:.2f}.".format(label = spectral_label, ami = ami_louvain), "yellow")
	