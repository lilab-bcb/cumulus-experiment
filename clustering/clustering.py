import numpy as np
import scCloud
from scipy.sparse.csgraph import connected_components
from scipy.sparse import issparse
from scCloud.tools.diffusion_map import calculate_affinity_matrix

def run_spectral(data, rep_key, K = 100, n_jobs = 1, random_state = 0, full_speed = False):

	indices, distances = get_kNN(data, rep_key, K, n_jobs = n_jobs, random_state = random_state, full_speed = full_speed)
	W = calculate_affinity_matrix(indices, distances)

	assert issparse(W)

	Phi_pt, Lambda = calculate_diffusion_map(W, n_dc = n_components, alpha = alpha, solver = solver, random_state = random_state)

	nc, labels = connected_components(W, directed = True, connection = 'strong')

	assert nc == 1

