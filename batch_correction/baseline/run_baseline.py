import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("scCloud cluster -p {cores} --run-umap ../../MantonBM_nonmix_tiny_10x.h5 baseline_result".format(cores = n_jobs)):
	sys.exit(1)