import scCloud
import os, sys
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_10x.h5"
data_dst = "MantonBM_nonmix_clustering"

def calc_and_plot():
	if os.system("scCloud cluster -p {jobs} --run-louvain --run-approximated-louvain --run-leiden --run-approximated-leiden --run-tsne {src} {outname}".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis tsne --attributes louvain_labels,Individual {name}.h5ad {name}.louvain_labels.tsne.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis tsne --attributes approx_louvain_labels,Individual {name}.h5ad {name}.approx_louvain_labels.tsne.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis tsne --attributes leiden_labels,Individual {name}.h5ad {name}.leiden_labels.tsne.pdf".format(name = data_dst)):
		sys.exit(1)

	if os.system("scCloud plot scatter --basis tsne --attributes approx_leiden_labels,Individual {name}.h5ad {name}.approx_leiden_labels.tsne.pdf".format(name = data_dst)):
		sys.exit(1)

def measure_results():
	adata = scCloud.tools.read_input(data_dst + '.h5ad', mode = 'a')

	cprint("Calculating AMI...", "green")
	ami_louvain = adjusted_mutual_info_score(adata.obs['louvain_labels'], adata.obs['approx_louvain_labels'], average_method = 'arithmetic')
	cprint("AMI between louvain_labels and approx_louvain_labels: {:.2f}.".format(ami_louvain), "yellow")
	ami_leiden = adjusted_mutual_info_score(adata.obs['leiden_labels'], adata.obs['approx_leiden_labels'], average_method = 'arithmetic')
	cprint("AMI between leiden_labels and approx_leiden_labels: {:.2f}.".format(ami_leiden), "yellow")

if __name__ == '__main__':
	if (data_dst + '.h5ad') not in os.listdir('.'):
		calc_and_plot()

	measure_results()
	