import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("sccloud cluster -p {cores} --umap ../../MantonBM_nonmix_tiny.h5sc baseline_result".format(cores = n_jobs)):
	sys.exit(1)

if os.system("sccloud de_analysis -p 8 --labels spectral_leiden_labels --t baseline_result.h5ad baseline_result.de.xlsx"):
	sys.exit(1)

if os.system("sccloud annotate_cluster baseline_result.h5ad baseline_result.anno.txt"):
	sys.exit(1)