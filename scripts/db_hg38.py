import pathlib as path
import src.scAnnotate.gff

if __name__ == '__main__':
    gff_path = path.Path(snakemake.input.gff)
    db_path = path.Path(snakemake.output.db)
    src.scAnnotate.gff.build_feature_db(gff_path, db_path)

