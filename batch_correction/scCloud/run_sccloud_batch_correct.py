import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("sccloud cluster -p {cores} --correct-batch-effect --umap ../../MantonBM_nonmix_tiny.h5sc tiny_sccloud_corrected".format(cores = n_jobs)):
	sys.exit(1)