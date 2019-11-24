import os, sys
import pegasus as pg
from termcolor import cprint

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"

basis_list = ['fitsne', 'tsne', 'net_tsne', 'umap', 'net_umap', 'fle', 'net_fle']

fig_index = {
	'fitsne'   : 'Figure_S6A_right',
	'tsne'     : 'Figure_S6A_left',
	'net_tsne' : 'Figure_S6A_center',
	'umap'     : 'Figure_2D_left',
	'net_umap' : 'Figure_2D_right',
	'fle'      : 'Figure_S6B_left',
	'net_fle'  : 'Figure_S6B_right'
}

palettes_str = "#c5b0d5,#ff7f0e,#8c564b,#ff9896,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#9edae5,#aec7e8,#ffbb78,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000"

def calc_and_plot(in_file):

	for basis in basis_list:
		if os.system('pegasus plot scatter --basis {basis} --attributes anno_louvain --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/{dst}.pdf'.format(basis = basis, palettes = palettes_str, src = in_file, dst = fig_index[basis])):
			sys.exit(1)


def measure_algorithms(adata):
	cprint("Calculating kSIM measurements...", "green")

	for basis in basis_list:
		ksim, ac_rate = pg.calc_kSIM(adata, attr = 'louvain_labels', rep = basis, n_jobs = n_cores)
		cprint("For {basis}, kSIM is {ksim:.4f}, accept rate is {ac_rate:.4f}.".format(basis = basis, ksim = ksim, ac_rate = ac_rate), "yellow")

def net_umap_phase_plots(adata):
	cprint("Generate phase plots for Net UMAP...", "green")

	# Initial plot
	bdata = adata[adata.obs['ds_selected'], :]
	bdata.obsm['X_net_umap_init'] = bdata.uns['X_net_umap_small']
	pg.write_output(bdata, "net_umap_init")

	if os.system('pegasus plot scatter --basis net_umap_init --attributes louvain_labels --set-palettes "{palettes}" net_umap_init.h5ad /output/Figure_2C_left.pdf'.format(palettes = palettes_str)):
		sys.exit(1)

	# Predict plot
	if os.system('pegasus plot scatter --basis net_umap_pred --attributes louvain_labels --set-palettes "{palettes}" {src}.h5ad /output/Figure_2C_center.pdf'.format(src = data_src, palettes = palettes_str)):
		sys.exit(1)


if __name__ == '__main__':
	calc_and_plot(data_src)

	adata = pg.read_input(data_src + '.h5ad')

	measure_algorithms(adata)
	net_umap_phase_plots(adata)