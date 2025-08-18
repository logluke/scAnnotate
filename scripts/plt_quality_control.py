import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == "__main__":
    data = sc.read_h5ad(snakemake.input.h5ad)

    # Plot qc metrics
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=[6.4 * 2 / 3, 4.8 * 3 / 2], tight_layout=True)
    ax2.sharex(ax1)
    fig.suptitle('Cell QC metrics')
    X = data.obs[['log1p_total_counts', 'log1p_total_genes', 'pct_mt_counts', 'time', 'cell_qc']]
    sns.histplot(X, x="log1p_total_genes", hue="cell_qc", ax=ax1, bins=25, multiple="stack")
    sns.histplot(X, x="log1p_total_counts", hue="cell_qc", ax=ax2, bins=25, multiple="stack")
    sns.histplot(X, x="pct_mt_counts", hue="cell_qc", ax=ax3, bins=25, multiple="stack")
    plt.savefig(snakemake.output.qc_marginal_png)

    h = sns.jointplot(X, x="log1p_total_counts", y="log1p_total_genes", kind="hist", hue="cell_qc", alpha=.75)
    h.ax_joint.set_xlabel('log$_{10}$ counts per cell')
    h.ax_joint.set_ylabel('log$_{10}$ genes per cell')

    N = 200
    X = np.linspace(2.0, 6.0, N)
    Y = np.linspace(1.5, 4.0, N)
    X, Y = np.meshgrid(X, Y)
    pos = np.dstack((X, Y))

    plt.savefig(snakemake.output.qc_joint_png)
    plt.show()
