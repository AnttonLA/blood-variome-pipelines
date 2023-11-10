import os
import sys
import polars as pl
from variant_list_to_fabian_input_vcf import variant_list_to_fabian_input_vcf, create_map_file


def test_variant_list_to_fabian_input_vcf():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    snplist = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv")
    output_dir = os.path.join(tests_dir, "tmp/")
    variant_list_to_fabian_input_vcf(snplist, output_dir)

    # Assert that at least one output file exist
    assert os.path.isfile(os.path.join(tests_dir, "tmp/FABIAN_INPUT_1.vcf"))

    # Assert that the output file has the correct format
    df = pl.read_csv(os.path.join(tests_dir, "tmp/FABIAN_INPUT_1.vcf"), separator='\t', has_header=True)
    assert df.columns == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA00001"]

    # Delete files once test is done
    os.remove(os.path.join(tests_dir, "tmp/FABIAN_INPUT_1.vcf"))


def test_create_map_file():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    snplist = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv")
    map_file = os.path.join(tests_dir, "tmp/dummy_variant_list.tsv.map")
    reference = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv.map")
    df = pl.read_csv(snplist, separator='\t', has_header=True)
    create_map_file(df, map_file)

    # Check that the map file exists
    assert os.path.isfile(map_file)

    # Check that the map file is the same as the reference file
    df = pl.read_csv(map_file, separator='\t', has_header=True)
    df_ref = pl.read_csv(reference, separator='\t', has_header=True)

    assert not df.is_empty()
    assert df.frame_equal(df_ref)
