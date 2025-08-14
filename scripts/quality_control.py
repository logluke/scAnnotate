import scanpy as sc
import src.scAnnotate.qc

if __name__ == "__main__":
    data = sc.read_h5ad(snakemake.input.h5ad)
    # Add standard annotations (log)total_genes, (log)total_counts, (log)pct_mt_counts
    src.scAnnotate.qc.add_annotations(data)

    # Cluster cells into "pass" and "fail" using GMM based on previous annotations.
    src.scAnnotate.qc.cell_quality_control(data)

    data.write_h5ad(snakemake.output.h5ad_qc, compression="gzip")
