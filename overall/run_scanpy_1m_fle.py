import time
import scanpy as sc
from datetime import timedelta

sc.settings.verbosity = 4
rand_seed = 0

adata = sc.read_h5ad("scanpy_neuron_result_no_fle.h5ad")
print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, iterations = 5000)
end_fle = time.time()
logstr_fle = "Time spent for FLE with {iters} iterations: {time}.".format(iters = 500, time = timedelta(seconds = end_fle - start_fle))
print(logstr_fle)