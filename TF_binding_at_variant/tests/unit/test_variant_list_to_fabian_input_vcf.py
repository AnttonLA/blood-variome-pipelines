import os
import sys
import polars as pl
from variant_list_to_fabian_input_vcf import variant_list_to_fabian_input_vcf


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
