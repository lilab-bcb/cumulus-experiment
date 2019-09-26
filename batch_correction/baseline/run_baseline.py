import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("sccloud cluster -p {cores} --umap /projects/benchmark/MantonBM/MantonBM_nonmix_tiny.h5sc baseline_result".format(cores = n_jobs)):
	sys.exit(1)