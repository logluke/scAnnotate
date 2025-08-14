import pysam
import numpy as np
import pandas as pd
import anndata as ad
import itertools as it
import pathlib as path
import dataclasses as dc
import multiprocessing as mp

from tqdm import tqdm
from scipy import sparse


import src.scAnnotate.gff


def is_high_quality_flag(flag: int) -> bool:
    """ Check if the flag indicates a high-quality read. """
    unmapped = 0b0000_0000_0000_0100  # 4
    secondary = 0b0000_0001_0000_0000  # 256
    qcfail = 0b0000_0010_0000_0000  # 512
    duplicate = 0b0000_0100_0000_0000  # 1024

    # Operator |: Sets each bit to 1 if one of two bits is 1
    bad_mask = (unmapped | secondary | qcfail | duplicate)

    # Operator &: Sets each bit to 1 if both bits are 1
    return (flag & bad_mask) == 0


def init_worker(bam_path: str,
                gene2col: dict[str, int],
                cell2row: dict[str, int],
                include_introns: bool):
    global shared_bam_reader, shared_gene2col, shared_cell2row, shared_use_introns

    shared_bam_reader = pysam.AlignmentFile(bam_path, "rb")
    shared_gene2col = gene2col
    shared_cell2row = cell2row
    shared_use_introns = include_introns



def worker(gene) -> tuple[int, np.array]:
    bam = shared_bam_reader
    g2c = shared_gene2col
    c2r = shared_cell2row
    include_introns = shared_use_introns
    n_cells = len(c2r)
    gene_id = gene.id
    col = g2c[gene_id]

    vec = np.zeros(n_cells, dtype=int)

    if include_introns:
        regions = [(gene.start, gene.end)]
    else:
        regions = gene.ivs

    for start, end in regions:
        for aln in bam.fetch(gene.chrom, start, end):
            if not is_high_quality_flag(aln.flag):
                continue
            if aln.mapping_quality < 30:
                continue
            if gene.strand== "+" and aln.is_reverse:
                continue
            if gene.strand== "-" and aln.is_forward:
                continue
            try:
                cb = aln.get_tag("CB")
            except KeyError:
                continue
            idx = c2r.get(cb)
            if idx is not None:
                vec[idx] += 1
    return col, vec



def count_matrix(bam_path : path.Path,
                 db_path : path.Path,
                 cm_path : path.Path,
                 *,
                 chrom : str = "all",
                 include_introns : bool = True):

    # Generate feature (gene > transcript > exon) tree to count reads
    genes = src.scAnnotate.gff.load_genes_with_exons(db_path, chrom)

    # Init AnnData Count Matrix to store coverage
    # Gene id as row label
    gene_ids = [gene.id for gene in genes]
    # Cell barcode as column label
    cell_ids = [f"{i:02d}_{j:02d}_{k:02d}" for i, j, k in it.product(range(1, 49), range(1, 97), range(1, 97))]

    # asign idx to id for quick access later
    gene2col = {gene_id: idx for idx, gene_id in enumerate(gene_ids)}
    cell2row = {cell_id: idx for idx, cell_id in enumerate(cell_ids)}
    n_genes = len(gene_ids)
    n_cells = len(cell_ids)
    global_counts = sparse.lil_matrix((n_cells, n_genes), dtype=int)

    # Distribute read processing for each gene to one CPU.
    with mp.Pool(
        processes=mp.cpu_count(),
        initializer=init_worker,
        initargs=(bam_path, gene2col, cell2row, include_introns),
    ) as pool:
        # map over all gene_ids, with a progress bar
        for col, vec in tqdm(
            pool.imap_unordered(worker, genes), total=len(genes)
        ):
            # fold each gene-vector into the global sparse matrix
            global_counts[:, col] = sparse.csr_matrix(vec).T

    # Convert the sparse matrix to a CSR format and create an AnnData object for storage
    counts = global_counts.tocsr()

    # AnnData reads genes in columns and cells in rows
    adata = ad.AnnData(
        X=counts, obs=pd.DataFrame(index=cell_ids), var=pd.DataFrame(index=gene_ids)
    )

    # Adding gene annotations
    var_annotation = pd.DataFrame.from_records((dc.asdict(g) for g in genes), index="id")
    if not include_introns:
        var_annotation["exon_start"] = var_annotation["ivs"].apply(lambda gene: [iv[0] for iv in gene])
        var_annotation['exon_start'] = var_annotation["exon_start"].apply(lambda x: ', '.join([str(i) for i in x]))
        var_annotation["exon_end"] = var_annotation["ivs"].apply(lambda gene: [iv[1]for iv in gene])
        var_annotation['exon_end'] = var_annotation["exon_end"].apply(lambda x: ', '.join([str(i) for i in x]))
    var_annotation.drop(columns=["exons"], inplace=True)
    adata.var = adata.var.join(var_annotation)

    # Adding cell metadata
    obs_annotation = pd.DataFrame({"time": np.repeat(["T0", "T1", "T2"], 16 * 96 * 96)}, index=adata.obs_names)
    adata.obs = adata.obs.join(obs_annotation)

    adata.write_h5ad(cm_path, compression="gzip")