Parse aligned single cell RNA-seq data (.bam) to count matrix (.h5ad).

# Usage

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install snakemake
```

# Run

```bash
snakemake --cores 1
```