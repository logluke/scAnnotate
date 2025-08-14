rule all:
    input:
        "data/gff/Homo_sapiens.GRCh38.114.chr.gff3.gz",
        "data/gff/Homo_sapiens.GRCh38.114.chr.db",
        "data/bam/example.chr21.bam",
        "data/h5ad/example.chr21.h5ad",
        "data/h5ad/example.chr21.qc.h5ad"


rule download_hg38:
    output:
        gff = "data/gff/Homo_sapiens.GRCh38.114.chr.gff3.gz"
    script:
        "scripts/download_hg38.py"


rule database_hg38:
    input:
        gff = "data/gff/Homo_sapiens.GRCh38.114.chr.gff3.gz"
    output:
        db = "data/gff/Homo_sapiens.GRCh38.114.chr.db"
    script:
        "scripts/database_hg38.py"


rule count_matrix:
    input:
        bam = "data/bam/example.chr21.bam",
        db = "data/gff/Homo_sapiens.GRCh38.114.chr.db"
    output:
        h5ad = temp("data/h5ad/example.chr21.h5ad")
    params:
        chrom = "21",
        include_introns = True
    script:
        "scripts/count_matrix.py"

rule quality_control:
    input:
        h5ad = "data/h5ad/example.chr21.h5ad"
    output:
        h5ad_qc = protected("data/h5ad/example.chr21.qc.h5ad")
    script:
        "scripts/quality_control.py"



