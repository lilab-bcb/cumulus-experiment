import os, sys
import pegasus as pg

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"
out_name = "MantonBM_nonmix_pegasus"

basis_list = ['fitsne', 'tsne', 'net_tsne', 'umap', 'net_umap', 'fle', 'net_fle']

def calc_and_plot(in_file, out_file):
	cprint("Calculation with {jobs} cores.".format(jobs = n_cores), "green")

	for basis in basis_list:
		if os.system("sccloud plot scatter --basis {basis} --attributes louvain_labels {src}.h5ad {dst}.louvain.{basis}.pdf".format(basis = basis, src = in_file, dst = out_file)):
			sys.exit(1)


def measure_algorithms(adata):
	cprint("Calculating kSIM measurements...", "green")

	for basis in basis_list:
		ksim, ac_rate = pg.calc_kSIM(adata, attr = 'louvain_labels', rep = basis, n_jobs = n_cores)
		print("For {basis}, kSIM is {ksim:.4f}, accept rate is {ac_rate:.4f}.".format(basis = basis, ksim = ksim, ac_rate = ac_rate))

def net_umap_phase_plots(adata, out_file):
	cprint("Generate phase plots for Net UMAP...", "green")

	# Initial plot
	bdata = adata[adata.obs['ds_selected'], :]
	bdata.obsm['X_net_umap_init'] = bdata.uns['X_net_umap_small']
	pg.write_output(bdata, "net_umap_init")

	if os.system("pegasus plot scatter --basis net_umap_init --attributes louvain_labels net_umap_init.h5ad {name}.net_umap_init.pdf".format(name = out_file)):
		sys.exit(1)

	# Predict plot
	if os.system("pegasus plot scatter --basis net_umap_pred --attributes louvain_labels {src}.h5ad {dst}.net_umap_pred.pdf".format(src = data_src, dst = out_file)):
		sys.exit(1)


if __name__ == '__main__':
	calc_and_plot(data_src, out_name)

	adata = pg.read_input(data_src + '.h5ad')

	measure_algorithms(adata)
	net_umap_phase_plots(adata, out_name)