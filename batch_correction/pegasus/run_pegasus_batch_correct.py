import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("pegasus cluster -p {cores} --correct-batch-effect --umap /data/MantonBM_nonmix_tiny.h5sc tiny_sccloud_corrected > pegasus_correct.log".format(cores = n_jobs)):
	sys.exit(1)