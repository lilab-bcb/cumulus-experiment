import sys, os
import numpy as np
import scCloud
from scCloud.tools.nearest_neighbors import calculate_nearest_neighbors
from distance_metric import get_kDistances
from termcolor import cprint

scanpy_methods = ["combat", "mnn"]

def run_exp(method, dataset):
	if method == 'scCloud':
		attr = "Channel"
	else:
		attr = "channel"

	# For dataset not batch corrected.
	if method == 'scCloud':
		if dataset == "tiny":
			src_file = './scCloud/MantonBM_nonmix_tiny_sccloud.h5ad'
			fn_prefix = 'MantonBM_nonmix_tiny_sccloud'
		else:
			src_file = 'MantonBM_nonmix.h5ad'
			fn_prefix = 'MantonBM_nonmix'
	else:
		if dataset == "tiny":
			src_file = './scanpy/MantonBM_nonmix_tiny_scanpy.h5ad'
			fn_prefix = 'MantonBM_nonmix_tiny_scanpy'
		else:
			src_file = 'MantonBM_nonmix_scanpy.h5ad'
			fn_prefix = 'MantonBM_nonmix_scanpy'

	adata = scCloud.tools.read_input(src_file, mode = 'a')
	if method in scanpy_methods:
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]

	kbet1, pvalue1, ac_rate1 = get_kDistances(adata, fn_prefix, attr)

	# For dataset batch corrected.
	if method == 'scCloud':
		if dataset == "tiny":
			src2_file = './scCloud/MantonBM_nonmix_tiny_sccloud_corrected.h5ad'
			fn2_prefix = 'MantonBM_nonmix_tiny_sccloud_corrected'
		else:
			src2_file = 'MantonBM_nonmix_corrected.h5ad'
			fn2_prefix = 'MantonBM_nonmix_corrected'
	elif method == 'combat':
		if dataset == "tiny":
			src2_file = './scanpy/combat/MantonBM_nonmix_tiny_scanpy_combat.h5ad'
			fn2_prefix = 'MantonBM_nonmix_tiny_scanpy_combat'
		else:
			src2_file = 'MantonBM_nonmix_scanpy_combat.h5ad'
			fn2_prefix = 'MantonBM_nonmix_scanpy_combat'
	else:
		if dataset == "tiny":
			src2_file = './scanpy/mnn/MantonBM_nonmix_tiny_scanpy_mnn_processed.h5ad'
			fn2_prefix = 'MantonBM_nonmix_tiny_scanpy_mnn'
		else:
			src2_file = 'MantonBM_nonmix_scanpy_mnn.h5ad'
			fn2_prefix = 'MantonBM_nonmix_scanpy_mnn'
		
	adata = scCloud.tools.read_input(src2_file, mode = 'a')	
	if method in scanpy_methods:
		adata.uns['knn_indices'] = adata.uns['neighbors']['indices'][:, 1:]	

	kbet2, pvalue2, ac_rate2 = get_kDistances(adata, fn_prefix, attr)

	cprint("kBET is improved by {:.2%}.".format(np.abs(kbet1 - kbet2) / kbet1))
	cprint("pvalue of kBET is improved by {:.2%}.".format(np.abs(pvalue1 - pvalue2) / pvalue1))
	cprint("Accept rate is improved by {:.2%}".format(np.abs(ac_rate1 - ac_rate2) / ac_rate1))


if __name__ == '__main__':
	method = sys.argv[1]
	dataset = sys.argv[2]
	assert method in (["scCloud"] + scanpy_methods) and dataset in ["tiny", "complete"]
	run_exp(method, dataset)