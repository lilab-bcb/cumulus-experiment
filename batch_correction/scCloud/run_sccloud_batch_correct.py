import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("scCloud cluster -p {cores} --correct-batch-effect --run-leiden --run-approximated-leiden --run-fitsne --run-umap ../../MantonBM_nonmix_tiny_10x.h5 tiny_sccloud_corrected".format(cores = n_jobs)):
	sys.exit(1)

if os.system("scCloud de_analysis -p {cores} --labels leiden_labels tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.de.xlsx".format(cores = n_jobs)):
	sys.exit(1)

if os.system("scCloud annotate_cluster tiny_sccloud_corrected.h5ad tiny_sccloud_corrected.anno.txt"):
	sys.exit(1)