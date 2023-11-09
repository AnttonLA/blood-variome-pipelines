import os
import polars as pl
from remap_lookup_for_full_snplist import remap_lookup_for_full_snplist


def test_remap_lookup_for_full_snplist():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    snplist = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv")
    remap_file = os.path.join(tests_dir, "test_data/dummy_remap_file.bed.gz")
    tmp_dir = os.path.join(tests_dir, "tmp/")

    # Make a list of the expected output files
    df = pl.read_csv(snplist, separator="\t", has_header=True)
    df = df.with_columns([pl.format("{}:{}", pl.col("Chrom").str.replace("chr", ""), pl.col("Pos")).alias("chr_pos")])
    df = df.with_columns((pl.lit("remap_studies_") + pl.col("chr_pos") + pl.lit(".txt")).alias("filename"))
    expected_output_files = df.get_column("filename").to_list()

    remap_lookup_for_full_snplist(snplist, remap_file, tmp_dir, tmp_dir)  # Output to tmp_dir

    # Assert that the output files exist
    for file in expected_output_files:
        assert os.path.isfile(os.path.join(tests_dir, "tmp/" + file))

    # Delete files once test is done
    for file in expected_output_files:
        os.remove(os.path.join(tests_dir, "tmp/" + file))
