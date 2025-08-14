rule all:
    input:
        "data/in/Homo_sapiens.GRCh38.114.chr.gff3.gz"

rule download_hg38:
    output:
        gff = "data/in/Homo_sapiens.GRCh38.114.chr.gff3.gz"
    script:
        "scripts/dl_hg38.py"