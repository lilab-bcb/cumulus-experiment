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

src_file = "/data/MantonBM_nonmix.h5sc"
pegasus_src = "../MantonBM_nonmix_pegasus"
outname = "MantonBM_pegasus_output"
n_cores = os.cpu_count()

palette_diffmap = "#c5b0d5,#ff7f0e,#8c564b,#ff9896,#1f77b4,#dbdb8d,#e377c2,#2ca02c,#9edae5,#aec7e8,#ffbb78,#98df8a,#d62728,#9467bd,#c49c94,#f7b6d2,#bcbd22,#17becf,#ad494a,#8c6d31,#000000"
anno_str = "anno_louvain:louvain_labels:1. CD4+ Naive T cells;2. CD14+ Monocytes;3. T helper cells;4. Naive B cells;5. Cytotoxic T cells;6. CD8+ Naive T cells;7. NK cells;8. Cytotoxic T cells;9. Erythrocytes;10. Memory B cells;11. CD14+ Monocytes;12. Pre B cells;13. HSCs;14. cDCs;15. CD16+ Monocytes;16. Pro B cells;17. pDCs;18. Plasma cells;19. Erythrocytes;20. Megakaryocytes;21. MSCs"

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
    ax.get_figure().savefig("/output/Figure_S4B.pdf")
    plt.close()

def gen_fig_s4a():
    ndc_list = [(15, 'left'), (50, 'center'), (100, 'right')]
    
    for ndc, position in ndc_list:
        if os.system("pegasus cluster -p {jobs} --correct-batch-effect --diffmap-ndc {num} --diffmap-maxt -1 --louvain --fle {src} MantonBM_nonmix_ndc_{num}_t_neg_one".format(jobs = n_cores, num = ndc, src = src_file)):
            sys.exit(1)

        if os.system('pegasus annotate_cluster --annotation "{anno}" MantonBM_nonmix_ndc_{num}_t_neg_one'.format(anno = anno_str, num = ndc)):
            sys.exit(1)

        if os.system('pegasus plot scatter --basis fle --attributes anno_louvain --wspace 1.2 --set-palettes "{palettes}" MantonBM_nonmix_ndc_{num}_t_neg_one.h5ad /output/Figure_S4A_{pos}.pdf'.format(palettes = palette_diffmap, pos = position)):
            sys.exit(1)


def gen_fig_s4b():
    rep = update_rep("pca")
    adata = pg.read_input("{}.h5ad".format(pegasus_src))
    t_arr, entropy_arr, knee_pt = get_entropy(W_from_rep(adata, rep))
    plot_curve(t_arr, entropy_arr, knee_pt)

def gen_fig_2b():
    if os.system("cp Figure_S4A_right.pdf /output/Figure_2B_left.pdf"):
        sys.exit(1)

    if os.system('pegasus plot scatter --basis fle --attributes anno_louvain --wspace 1.2 --set-palettes "{palettes}" {src}.h5ad /output/Figure_2B_right.pdf'.format(palettes = palette_diffmap, src = pegasus_src)):
        sys.exit(1)

if __name__ == '__main__':
    gen_fig_s4a()
    gen_fig_s4b()
    gen_fig_2b()
    