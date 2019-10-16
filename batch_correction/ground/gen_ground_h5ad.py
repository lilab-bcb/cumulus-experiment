import os, sys

n_cores = os.cpu_count()

if os.system("pegasus cluster -p {} --louvain --umap /data/MantonBM_nonmix_tiny.h5sc ground".format(n_cores)):
	sys.exit(1)

if os.system("pegasus de_analysis -p {} --labels louvain_labels --t ground.h5ad ground.de.xlsx".format(n_cores)):
	sys.exit(1)

if os.system("pegasus annotate_cluster ground.h5ad ground.anno.txt"):
	sys.exit(1)

if os.system("pegasus plot scatter --basis umap --attributes louvain_labels,Channel ground.h5ad ground.umap.pdf"):
	sys.exit(1)