import scCloud
import numpy as np
import pandas as pd
import os, sys
from termcolor import cprint

markers = np.array(['CD3D', 'CD3E', 'CD3G', 'TRAC', 'CD4', 'CD8A', 'CD8B', 'RTKN2', 'TIGIT', 'FOXP3', 'IL2RA', 'CCR7', 'SELL', 'IL7R', 'TCF7', 'CD27', 'CD19', 'MS4A1', 'CD79A', 'CD79B', 'MZB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'NCAM1', 'NKG7', 'KLRB1', 'KLRD1', 'KLRF1', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRC4', 'FCGR3A', 'FCGR3B', 'ITGAL', 'ITGAM', 'VCAN', 'FCN1', 'S100A8', 'S100A9', 'CD14', 'CSF3R', 'CSF1R', 'CX3CR1', 'ITGAX', 'HLA-DPB1', 'HLA-DPA1', 'HLA-DQA1', 'CD1C', 'CD1E', 'FCER1A', 'CLEC10A', 'FCGR2B', 'THBD', 'CLEC9A', 'CADM1', 'XCR1', 'IL3RA', 'GZMB', 'JCHAIN', 'IRF7', 'TCF4', 'LILRA4', 'CLEC4C', 'CD38', 'XBP1', 'SLAMF7', 'CD34', 'KIT', 'CD59', 'THY1', 'SOX4', 'GYPA', 'HBB', 'HBA1', 'TFRC', 'ITGA4', 'ANK1', 'ICAM4', 'SLC4A1', 'ACKR1', 'BCAM', 'PF4', 'PPBP', 'GP5', 'SLAMF1', 'MPL', 'CXCR4', 'ITGA2B', 'CEACAM8', 'FUT4', 'MPO'])



if __name__ == '__main__':
	f_list = [f for f in os.listdir('.') if f in ["MantonBM_nonmix_seurat_corrected.h5ad", "MantonBM_nonmix_new_corrected.h5ad"]]
	if len(f_list) != 2:
		cprint("Computing highly variable genes using Seurat-style method...", "green")
		if os.system("scCloud cluster -p 8 --correct-batch-effect --select-hvg-flavor Seurat --run-leiden --run-tsne --run-fitsne --run-umap --run-fle ../MantonBM_nonmix_10x.h5 MantonBM_nonmix_seurat_corrected"):
			sys.exit(1)

		cprint("Computing highly variable genes using scCloud new method...", "green")
		if os.system("scCloud cluster -p 8 --plot-hvg --correct-batch-effect --run-leiden --run-tsne --run-fitsne --run-umap --run-fle ../MantonBM_nonmix_10x.h5 MantonBM_nonmix_new_corrected"):
			sys.exit(1)

