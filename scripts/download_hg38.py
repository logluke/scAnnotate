import pathlib as path
import src.scAnnotate.gff

if __name__ == "__main__":
    gff_path = path.Path(snakemake.output.gff)
    # Download the Human reference genome GFF file from Ensembl
    src.scAnnotate.gff.download_gff(gff_path)