import scCloud
import os, sys
import numpy as np
from termcolor import cprint
from scCloud.tools.visualization import calc_umap

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_10x.h5"
data_dst = "MantonBM_nonmix_origin_vs_net"

basis_list = ['fitsne', 'tsne', 'net_tsne', 'umap', 'net_umap', 'fle', 'net_fle']

def calc_and_plot():
	cprint("Calculation with {jobs} cores.".format(jobs = n_cores), "green")
	
	if os.system("scCloud cluster -p {jobs} --correct-batch-effect --run-approximated-leiden --run-fitsne --run-tsne --run-net-tsne --run-umap --run-net-umap --run-fle --run-net-fle {src} {outname}".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

	for basis in basis_list:
		if os.system("scCloud plot scatter --basis {basis} --attributes approx_leiden_labels {name}.h5ad {name}.approx_leiden_labels.{basis}.pdf".format(basis = basis, name = data_dst)):
			sys.exit(1)

	if os.system("scCloud de_analysis -p {jobs} --labels approx_leiden_labels {name}.h5ad {name}.de.xlsx".format(jobs = n_cores, name = data_dst)):
		sys.exit(1)

	if os.system("scCloud annotate_cluster {name}.h5ad {name}.anno.txt".format(name = data_dst)):
		sys.exit(1)


def measure_algorithms(adata):
	cprint("Calculating kSIM measurements...", "green")

	for basis in basis_list:
		ksim, ac_rate = scCloud.tools.calc_kSIM(adata, 'approx_leiden_labels', rep_key = 'X_' + basis, n_jobs = n_cores)
		cprint("For {basis}, kSIM is {ksim:.4f}, accept rate is {ac_rate:.4f}.".format(basis = basis, ksim = ksim, ac_rate = ac_rate), "yellow")

def net_umap_phase_plots(adata):
	cprint("Generate phase plots for Net UMAP...", "green")

	# Initial plot
	bdata = adata[adata.obs['ds_selected'], :]
	bdata.obsm['X_net_umap_init'] = bdata.uns['X_net_umap_small']
	scCloud.tools.write_output(bdata, "net_umap_init")

	if os.system("scCloud plot scatter --basis net_umap_init --attributes approx_leiden_labels net_umap_init.h5ad net_umap_init.umap.pdf"):
		sys.exit(1)

	# Predict plot
	if os.system("scCloud plot scatter --basis net_umap_pred --attributes approx_leiden_labels {name}.h5ad {name}.net_umap_pred.pdf".format(name = data_dst)):
		sys.exit(1)


if __name__ == '__main__':
	if (data_dst + '.h5ad') not in os.listdir('.'):
		calc_and_plot()

	adata = scCloud.tools.read_input(data_dst + '.h5ad', mode = 'a')

	#measure_algorithms(adata)
	net_umap_phase_plots(adata)