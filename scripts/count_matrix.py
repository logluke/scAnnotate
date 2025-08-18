import src.scAnnotate.cm

if __name__ == "__main__":
    # Generate count matrix
    src.scAnnotate.cm.count_matrix(
        bam_path=snakemake.input.bam,
        db_path=snakemake.input.db,
        cm_path=snakemake.output.h5ad,
        chrom=snakemake.params.chrom,
        include_introns=snakemake.params.include_introns,
    )
