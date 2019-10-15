import pegasus as pg
import os, sys
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"
out_name = "MantonBM_nonmix_pegasus"

label_list = ['louvain', 'spectral_louvain', 'leiden', 'spectral_leiden']


def gen_plots(in_file, out_file):

	for label in label_list:
		if os.system("pegasus plot scatter --basis fitsne --attributes {label}_labels {src}.h5ad {dst}.{label}.fitsne.pdf".format(label = label, src = in_file, dst = out_file)):
			sys.exit(1)


def measure_results():
	adata = pg.read_input(data_src + '.h5ad')

	cprint("Calculating AMI...", "green")
	ami_louvain = adjusted_mutual_info_score(adata.obs['louvain_labels'], adata.obs['spectral_louvain_labels'], average_method = 'arithmetic')
	cprint("AMI between louvain_labels and spectral_louvain_labels: {:.2f}.".format(ami_louvain), "yellow")
	ami_leiden = adjusted_mutual_info_score(adata.obs['leiden_labels'], adata.obs['spectral_leiden_labels'], average_method = 'arithmetic')
	cprint("AMI between leiden_labels and spectral_leiden_labels: {:.2f}.".format(ami_leiden), "yellow")

if __name__ == '__main__':

	gen_plots(data_src, out_name)
	measure_results()
	