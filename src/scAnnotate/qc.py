import numpy as np
import anndata as ad
from sklearn import mixture
import scanpy as sc

def add_annotations(data: ad.AnnData):
    # Calculate cell quality controll metrics
    n_genes_per_cell = (data.X > 1).sum(axis=1).A1
    data.obs['total_genes'] = n_genes_per_cell
    data.obs['log1p_total_genes'] = np.log10(n_genes_per_cell + 1)

    total_count_per_cell = data.X.sum(axis=1).A1
    data.obs['total_counts'] = total_count_per_cell
    data.obs['log1p_total_counts'] = np.log10(total_count_per_cell + 1)

    # Calculate percentage of mitochondrial gene coverage per cell
    mt_mask = data.var["chrom"] == "MT"
    mt_data = data[:, mt_mask]
    mt_count_per_cell = mt_data.X.sum(axis=1).A1
    data.obs['mt_counts'] = mt_count_per_cell + 1
    data.obs['pct_mt_counts'] = np.where(total_count_per_cell > 0, mt_count_per_cell / total_count_per_cell * 100, 0)
    data.obs['log_pct_mt_counts'] = np.log10(data.obs['pct_mt_counts'])
    # np.where(condition, x, y) picks x[i] where condition[i] is True, else y[i]

def cell_quality_control(data: ad.AnnData):
    X = data.obs[['log1p_total_counts', 'log1p_total_genes', 'pct_mt_counts']]

    # Fit GMM
    gmm = mixture.GaussianMixture(n_components=2, covariance_type="full")
    gmm.fit(X.values)  # DataFrame is array‐like :contentReference[oaicite:0]{index=0}

    means = gmm.means_
    g0_mean_total_count = means[0][0]
    g1_mean_total_count = means[1][0]

    if g0_mean_total_count > g1_mean_total_count:
        map = {0: "pass", 1: "fail"}
    else:
        map = {0: "fail", 1: "pass"}
    # Annotate cell quality control
    predictions = gmm.predict(X.values)
    data.obs["cell_qc"] = predictions
    data.obs["cell_qc"] = data.obs["cell_qc"].map(map)
