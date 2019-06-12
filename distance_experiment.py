import sys, os
import numpy as np
import scCloud
from scCloud.tools.nearest_neighbors import calculate_nearest_neighbors
from distance_metric import get_kDistances
from termcolor import cprint

def run_exp(package, attr, n_jobs = 8):
	print("Consider attribute " + attr + ":")

	# For dataset not batch corrected.
	if package == 'scCloud':
		adata = scCloud.tools.read_input('MantonBM_nonmix.h5ad', mode = 'a')
		fn_prefix = 'MantonBM_nonmix'
	else:
		adata = scCloud.tools.read_input('MantonBM_nonmix_scanpy.h5ad', mode = 'a')
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]
		fn_prefix = 'MantonBM_nonmix_scanpy'

	kbet1, pvalue1, kbjsd1 = get_kDistances(adata, fn_prefix, attr)

	# For dataset batch corrected.
	if package == 'scCloud':
		adata = scCloud.tools.read_input('MantonBM_nonmix_corrected.h5ad', mode = 'a')
		fn_prefix = 'MantonBM_nonmix_corrected'
	else:
		adata = scCloud.tools.read_input('MantonBM_nonmix_scanpy_combat.h5ad', mode = 'a')
		#indices, _ = calculate_nearest_neighbors(adata.obsm['X_pca'], num_threads = n_jobs)
		#adata.uns['knn_indices'] = indices
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]
		fn_prefix = 'MantonBM_nonmix_scanpy_combat'

	kbet2, pvalue2, kbjsd2 = get_kDistances(adata, fn_prefix, attr)

	cprint("kBET is improved by {:.2%}.".format(np.abs(kbet1 - kbet2) / kbet1))
	cprint("pvalue of kBET is improved by {:.2%}.".format(np.abs(pvalue1 - pvalue2) / pvalue1))
	cprint("kBJSD is improved by {:.2%}.".format(np.abs(kbjsd1 - kbjsd2) / kbjsd1))


if __name__ == '__main__':
	package = sys.argv[1]
	attr = sys.argv[2]
	assert package == 'scCloud' or package == 'scanpy'
	run_exp(package, attr)