import os
import sys
import gffutils
import urllib.request
import pathlib as path

from dataclasses import dataclass

Interval = tuple[int, int]

# @dataclass
# class Exon:
#     id: str
#     transcript_id: str
#     chrom: str
#     strand: str
#     start: int
#     end: int
#
#
# @dataclass
# class Transcript:
#     id: str
#     gene_id: str
#     chrom: str
#     strand: str
#     start: int
#     end: int
#     exons: list[Exon]


@dataclass
class Gene:
    id: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: list[Interval]


def download_gff(gff_path: path.Path | str):
    url = f"https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/{gff_path.name}"
    try:
        urllib.request.urlretrieve(url, gff_path)
        print(f"Download complete: {gff_path}")
    except Exception as e:
        print(f"Download failed: {e}")
        if os.path.exists(gff_path):
            os.remove(gff_path)
        sys.exit(1)


def build_feature_db(
    gff_path: path.Path,
    db_path: path.Path | str,
    *,
    force: bool = False,
    m_strat: str = "create_unique",
    keep_order: bool = True,
):
    gffutils.create_db(
        str(gff_path),
        dbfn=str(db_path),
        force=force,
        keep_order=keep_order,
        merge_strategy=m_strat,
    )


def load_feature_db(db_path: path.Path | str):
    try:
        db = gffutils.FeatureDB(str(db_path))
    except Exception as e:
        raise RuntimeError(f"Error loading GFF DB at {db_path}: {e}")
    return db


def merge_overlapping_intervals(intervals: list[Interval]) -> list[Interval]:
    if not intervals:
        return []

    # Validate
    for iv in intervals:
        if iv[0] > iv[1]:
            raise ValueError(f"Invalid interval with start > end: {iv}")

    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged: list[Interval] = []

    current_start, current_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        # If intervals overlap or touch
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return merged


def load_genes_with_exons(db_path: path.Path | str, chr_id: str) -> list[Gene]:
    db = load_feature_db(db_path)
    genes: list[Gene] = []

    for gene in db.features_of_type("gene"):
        if gene.chrom != chr_id and chr_id != "all":
            continue
        exon_ivs = [
            (exon.start - 1, exon.end)
            for exon in db.children(gene, featuretype="exon", level=2)
        ]
        m_exon_ivs = merge_overlapping_intervals(exon_ivs)
        gi = Gene(
            id=gene.id.replace("gene:", ""),
            chrom=gene.chrom,
            strand=gene.strand,
            start=gene.start-1,  # Convert to 0-based start
            end=gene.end,
            exons=m_exon_ivs,
        )
        genes.append(gi)
    return genes
