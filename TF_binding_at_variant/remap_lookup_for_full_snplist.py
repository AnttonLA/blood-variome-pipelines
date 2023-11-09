import os
import argparse
import polars as pl
from extract_remapdb_studies import extract_studies_for_single_snp

"""
Wrapper function to apply extract_studies_for_single_snp() to every SNP in a 'snplist' file.
"""


def remap_lookup_for_full_snplist(variant_list_file: str, remap_path: str, tmp_dir: str, output_dir: str) -> None:
    """
    This function will take a list of variants and produce a tabix query to the ReMap metadata file
    :param variant_list_file: file with the variants to look up
    :param remap_path: Path to ReMap metadata file. Requires tabix index file.
    :param tmp_dir: Directory to store temporary files (i.e. tabix "slices" of the remap metadata file for the queried
        position). If none is specified, temporary files will be stored in the current working directory.
    :param output_dir: Directory to store the output files
    :return:
    """
    # Read snplist
    df = pl.read_csv(variant_list_file, separator="\t", has_header=True,
                     dtypes={"ID": pl.Utf8, "Chrom": pl.Utf8, "Pos": pl.Int64, "OA": pl.Utf8, "EA": pl.Utf8})

    # Check that columns "ID", "Chrom", "Pos", "OA", "EA" exist
    if df.columns != ["ID", "Chrom", "Pos", "OA", "EA"]:
        raise ValueError("Input file does not have the correct columns. The columns should be: ID, Chrom, Pos, OA, EA")

    # Combine columns Chrom (without 'chr') and Pos with a ':' in between into a new column called chr_pos
    df = df.with_columns([pl.format("{}:{}", pl.col("Chrom").str.replace("chr", ""), pl.col("Pos")).alias("chr_pos")])

    # Run extract_studies_for_single_snp() for each chr_pos value. This will produce a file for each variant position.
    df = df.with_columns(pl.col("chr_pos")
                         .apply(lambda s: extract_studies_for_single_snp(s,  # s is the chr_pos value, <chr>:<pos>
                                                                         remap_path,
                                                                         tmp_dir,
                                                                         os.path.join(output_dir,
                                                                                      f"remap_studies_{str(s)}.txt"))))


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Run extract_studies_for_single_snp() for every SNP in a file")
    parser.add_argument("-s", "--snplist", required=True, help="Path to file with snplist")
    parser.add_argument("-r", "--remapdb", required=True, help="Path to ReMap database")
    parser.add_argument("-t", "--tmp_dir", required=True,
                        help="Path to temporary directory where intermediate files will be stored")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Path to output dir. If the snplist contains several SNPs, "
                             "the output will be as many files, "
                             "named remap_studies_<chr:pos>.txt")

    args = parser.parse_args()

    # Validate input
    if not os.path.isfile(args.snplist):
        raise ValueError("Input file does not exist")

    if not os.path.isdir(args.tmpdir):
        raise ValueError("Temporary directory does not exist")
