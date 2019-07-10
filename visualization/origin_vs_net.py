import scCloud
import os, sys
from termcolor import cprint

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_10x.h5"
data_dst = "MantonBM_nonmix_origin_vs_net"

basis_list = ['tsne', 'net_tsne', 'umap', 'net_umap', 'fle', 'net_fle']

def calc_and_plot():
	cprint("Calculation with {jobs} cores.".format(jobs = n_cores), "green")
	
	if os.system("scCloud cluster -p {jobs} --correct-batch-effect --run-approximated-leiden --run-tsne --run-net-tsne --run-umap --run-net-umap --run-fle --run-net-fle {src} {outname}".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

	for basis in basis_list:
		if os.system("scCloud plot scatter --basis {basis} --attributes approx_leiden_labels,Individual {name}.h5ad {name}approx_leiden_labels.{basis}.pdf".format(basis = basis, name = data_dst)):
			sys.exit(1)


def measure_algorithms():
	adata = scCloud.tools.read_input(data_dst + '.h5ad', mode = 'a')
	cprint("Calculating kSIM measurements...", "green")

	for basis in basis_list:
		ksim, ac_rate = scCloud.tools.calc_kSIM(adata, 'approx_leiden_labels', rep_key = 'X_' + basis, n_jobs = n_cores)
		cprint("For {basis}, kSIM is {ksim}, accept rate is {ac_rate}.".format(basis = basis, ksim = ksim, ac_rate = ac_rate), "yellow")

if __name__ == '__main__':
	if (data_dst + '.h5ad') not in os.listdir('.'):
		calc_and_plot()

	measure_algorithms()