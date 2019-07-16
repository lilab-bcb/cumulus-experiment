import scCloud
import os, sys, time
import numpy as np
from termcolor import cprint

n_cores = os.cpu_count()
output_dir = "sccloud_output"
out_name = "MantonBM_nonmix_sccloud_corrected"
data_src = "../MantonBM_nonmix_10x.h5"

if output_dir not in os.listdir('.'):
	if os.system("mkdir {}".format(output_dir)):
		sys.exit(1)

if (out_name + '.h5ad') not in os.listdir(output_dir):
	if os.system("scCloud cluster -p {jobs} --correct-batch-effect --diffmap-full-speed --run-approximated-leiden --run-fitsne --net-ds-full-speed --run-net-umap --run-net-fle {src_name} {out_name}".format(jobs = n_cores, src_name = data_src, out_name = output_dir + '/' + out_name)):
		sys.exit(1)

adata = scCloud.tools.read_input(output_dir + '/' + out_name + '.h5ad', mode = 'a')

if (out_name + '_hvg.h5ad') not in os.listdir(output_dir):
	adata_hvg = adata[:, adata.var['highly_variable_genes']].copy()
	scCloud.tools.write_output(adata_hvg, out_name + '_hvg')

if 'x_pca.txt' not in os.listdir(output_dir):
	X = adata.obsm['X_pca']
	np.savetxt(output_dir + '/x_pca.txt', X)
