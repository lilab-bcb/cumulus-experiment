import scCloud
import os, sys
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_10x.h5"
data_dst = "MantonBM_nonmix_clustering"

label_list = ['louvain_labels', 'approx_louvain_labels', 'leiden_labels', 'approx_leiden_labels']

def process_data():
	cprint("Calculation with {} cores.".format(n_cores), "green")

	if os.system("scCloud cluster -p {jobs} --correct-batch-effect --run-louvain --run-approximated-louvain --run-leiden --run-approximated-leiden --run-fitsne {src} {outname}".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

def gen_plots():

	for label in label_list:
		if os.system("scCloud plot scatter --basis fitsne --attributes {label} {name}.h5ad {name}.{label}.fitsne.pdf".format(label = label, name = data_dst)):
			sys.exit(1)


def annotate_clusters():

	for label in label_list:
		if os.system("scCloud de_analysis -p {jobs} --labels {label} {name}.h5ad {name}.{label}.de.xlsx".format(jobs = n_cores, label = label, name = data_dst)):
			sys.exit(1)

		if os.system("scCloud annotate_cluster {name}.h5ad {name}.{label}.anno.txt".format(name = data_dst, label = label)):
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
		process_data()

	gen_plots()
	annotate_clusters()

	measure_results()
	