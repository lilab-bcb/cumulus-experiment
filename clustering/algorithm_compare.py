import pegasus as pg
import os, sys
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"
out_name = "MantonBM_nonmix_pegasus"

label_list = ['louvain', 'spectral_louvain', 'leiden', 'spectral_leiden']

palettes_dict = {
	'louvain': "#c5b0d5,#ff7f0e,#8c564b,#ff9896,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#9edae5,#aec7e8,#ffbb78,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000",
	'leiden': "#c5b0d5,#ff7f0e,#ff9896,#8c564b,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#ffbb78,#9edae5,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#aec7e8,#8c6d31,#000000",
	'spectral_louvain': "#c5b0d5,#ff7f0e,#8c564b,#1f77b4,#ff9896,#dbdb8d,#e377c2,#2ca02c,#ffbb78,#aec7e8,#9edae5,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000",
	'spectral_leiden': "#c5b0d5,#ff7f0e,#8c564b,#ff9896,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#aec7e8,#ffbb78,#9edae5,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000"
}

def gen_plots(in_file):

	for label in label_list:
		if os.system('pegasus plot scatter --basis fitsne --attributes anno_{label} --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/Figure_S5.{label}.pdf'.format(label = label, src = in_file, palettes = palettes_dict[label])):
			sys.exit(1)


def measure_results():
	adata = pg.read_input(data_src + '.h5ad')

	cprint("Calculating AMI...", "green")
	ami_louvain = adjusted_mutual_info_score(adata.obs['louvain_labels'], adata.obs['spectral_louvain_labels'], average_method = 'arithmetic')
	cprint("AMI between louvain_labels and spectral_louvain_labels: {:.2f}.".format(ami_louvain), "yellow")
	ami_leiden = adjusted_mutual_info_score(adata.obs['leiden_labels'], adata.obs['spectral_leiden_labels'], average_method = 'arithmetic')
	cprint("AMI between leiden_labels and spectral_leiden_labels: {:.2f}.".format(ami_leiden), "yellow")

if __name__ == '__main__':

	gen_plots(data_src)
	measure_results()
	