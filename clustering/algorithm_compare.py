import sccloud as scc
import os, sys
from termcolor import cprint
from sklearn.metrics.cluster import adjusted_mutual_info_score

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix.h5sc"
data_dst = "MantonBM_nonmix_clustering"

label_list = ['louvain_labels', 'spectral_louvain_labels', 'leiden_labels', 'spectral_leiden_labels']

def process_data():
	cprint("Calculation with {} cores.".format(n_cores), "green")

	if os.system("sccloud cluster -p {jobs} --correct-batch-effect --diffmap --louvain --spectral-louvain --leiden --spectral-leiden --fitsne {src} {outname} > clustering.log".format(jobs = n_cores, src = data_src, outname = data_dst)):
		sys.exit(1)

def gen_plots():

	for label in label_list:
		if os.system("sccloud plot scatter --basis fitsne --attributes {label} {name}.h5ad {name}.{label}.fitsne.pdf".format(label = label, name = data_dst)):
			sys.exit(1)


def annotate_clusters():

	for label in label_list:
		if os.system("sccloud de_analysis -p {jobs} --labels {label} --t {name}.h5ad {name}.{label}.de.xlsx".format(jobs = n_cores, label = label, name = data_dst)):
			sys.exit(1)

		if os.system("sccloud annotate_cluster {name}.h5ad {name}.{label}.anno.txt".format(name = data_dst, label = label)):
			sys.exit(1)


def measure_results():
	adata = scc.read_input(data_dst + '.h5ad')

	cprint("Calculating AMI...", "green")
	ami_louvain = adjusted_mutual_info_score(adata.obs['louvain_labels'], adata.obs['spectral_louvain_labels'], average_method = 'arithmetic')
	cprint("AMI between louvain_labels and spectral_louvain_labels: {:.2f}.".format(ami_louvain), "yellow")
	ami_leiden = adjusted_mutual_info_score(adata.obs['leiden_labels'], adata.obs['spectral_leiden_labels'], average_method = 'arithmetic')
	cprint("AMI between leiden_labels and spectral_leiden_labels: {:.2f}.".format(ami_leiden), "yellow")

if __name__ == '__main__':
	if (data_dst + '.h5ad') not in os.listdir('.'):
		process_data()

	gen_plots()
	annotate_clusters()

	measure_results()
	