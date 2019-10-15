import os, sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pegasus as pg

from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.stats import entropy
from sklearn.utils.extmath import randomized_svd
from termcolor import cprint

from pegasus.tools import update_rep, W_from_rep
from pegasus.tools.diffusion_map import calculate_normalized_affinity, calc_von_neumann_entropy, find_knee_point

src_name = "../MantonBM_nonmix_pegasus"
outname = "MantonBM_pegasus_output"

def get_entropy(W: "csr_matrix", n_components: int = 100, solver: str = 'eigsh', random_state: int = 0, max_t: int = 5000):
    cprint("Calculate coordinates:", "yellow")

    assert issparse(W)

    nc, labels = connected_components(W, directed=True, connection="strong")

    assert nc == 1

    W_norm, diag, diag_half = calculate_normalized_affinity(W)

    if solver == "eigsh":
        np.random.seed(random_state)
        v0 = np.random.uniform(-1.0, 1.0, W_norm.shape[0])
        Lambda, U = eigsh(W_norm, k=n_components, v0=v0)
        Lambda = Lambda[::-1]
        U = U[:, ::-1]
    else:
        assert solver == "randomized"
        U, S, VT = randomized_svd(W_norm, n_components=n_components, random_state=random_state)
        signs = np.sign((U * VT.transpose()).sum(axis=0))  # get eigenvalue signs
        Lambda = signs * S  # get eigenvalues

    # remove the first eigen value and vector
    Lambda = Lambda[1:]
    U = U[:, 1:]
    Phi = U / diag_half[:, np.newaxis]

    # Find the knee point 
    x = np.array(range(1, max_t + 1), dtype = float)
    y = np.array([calc_von_neumann_entropy(Lambda, t) for t in x])
    t = x[find_knee_point(x, y)]

    return x.astype('int'), y, t

def plot_curve(t_arr, entropy_arr, knee_pt):
    cprint("Plotting:", "yellow")

    ax = sns.lineplot(x = t_arr, y = entropy_arr)
    ax = sns.scatterplot(x = np.array([knee_pt]), y = np.array([entropy_arr[int(knee_pt)]]), color = 'black')
    ax.text(knee_pt + 0.01, entropy_arr[int(knee_pt)] - 0.01, "Knee Point", horizontalalignment = 'left', size = 'medium', color = 'black')
    plt.xlabel("Time Stamp")
    plt.ylabel("von Neumann Entropy")
    ax.get_figure().savefig("Figure_S4B.pdf")
    plt.close()

if __name__ == '__main__':

    rep = update_rep("pca")
    adata = pg.read_input("{}.h5ad".format(src_name))
    t_arr, entropy_arr, knee_pt = get_entropy(W_from_rep(adata, rep))
    plot_curve(t_arr, entropy_arr, knee_pt)
