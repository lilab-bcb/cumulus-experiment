import scCloud
import os, sys
from termcolor import cprint

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_10x.h5"
data_dst = "MantonBM_nonmix_origin_vs_net"

def calc_and_plot():
	if os.system("scCloud cluster -p {jobs} --correct-batch-effect --run-approximated-leiden --run-tsne --run-net-tsne --run-umap --run-net-umap --run-fle --run-net-fle {src} {outname}".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis tsne --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.tsne.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis net_tsne --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.net_tsne.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis umap --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.umap.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis net_umap --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.net_umap.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis fle --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.fle.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis net_fle --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.net_fle.pdf".format(name = data_dst)):
		sys.exit(1)

def measure_algorithms():
	adata = scCloud.tools.read_input(data_dst + '.h5ad', mode = 'a')
	cprint("Calculating kSIM measurements...", "green")

	method_list = ['tsne', 'net_tsne', 'umap', 'net_umap', 'fle', 'net_fle']
	for method in method_list:
		ksim, ac_rate = scCloud.tools.calc_kSIM(adata, 'approx_leiden_labels', rep_key = 'X_' + method, n_jobs = n_cores)
		cprint("For {method}, kSIM is {ksim}, accept rate is {ac_rate}.".format(method = method, ksim = ksim, ac_rate = ac_rate), "yellow")

if __name__ == '__main__':
	if (data_dst + '.h5ad') not in os.listdir('.'):
		calc_and_plot()

	measure_algorithms()