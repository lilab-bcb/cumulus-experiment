import sys, os
import numpy as np
import scCloud
from scCloud.tools.nearest_neighbors import calculate_nearest_neighbors
from distance_metric import get_kDistances
from termcolor import cprint

def run_exp(method):
	if method == 'scCloud':
		attr = "Channel"
	else:
		attr = "channel"

	# For dataset not batch corrected.
	if method == 'scCloud':
		adata = scCloud.tools.read_input('MantonBM_nonmix.h5ad', mode = 'a')
		fn_prefix = 'MantonBM_nonmix'
	else:
		adata = scCloud.tools.read_input('./scanpy/MantonBM_nonmix_tiny_scanpy.h5ad', mode = 'a')
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]
		fn_prefix = 'MantonBM_nonmix_tiny_scanpy'

	kbet1, pvalue1, kbjsd1 = get_kDistances(adata, fn_prefix, attr)

	# For dataset batch corrected.
	if method == 'scCloud':
		adata = scCloud.tools.read_input('MantonBM_nonmix_corrected.h5ad', mode = 'a')
		fn_prefix = 'MantonBM_nonmix_corrected'
	else:
		if method == 'combat':
			adata = scCloud.tools.read_input('./scanpy/combat/MantonBM_nonmix_tiny_scanpy_combat.h5ad', mode = 'a')
			fn_prefix = 'MantonBM_nonmix_tiny_scanpy_combat'
		else:
			adata = scCloud.tools.read_input('./scanpy/bbknn/MantonBM_nonmix_tiny_scanpy_bbknn.h5ad', mode = 'a')
			fn_prefix = 'MantonBM_nonmix_tiny_scanpy_bbknn'
		
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]
		

	kbet2, pvalue2, kbjsd2 = get_kDistances(adata, fn_prefix, attr)

	cprint("kBET is improved by {:.2%}.".format(np.abs(kbet1 - kbet2) / kbet1))
	cprint("pvalue of kBET is improved by {:.2%}.".format(np.abs(pvalue1 - pvalue2) / pvalue1))
	cprint("kBJSD is improved by {:.2%}.".format(np.abs(kbjsd1 - kbjsd2) / kbjsd1))


if __name__ == '__main__':
	method = sys.argv[1]
	assert method in ["scCloud", "combat", "bbknn"]
	run_exp(method)