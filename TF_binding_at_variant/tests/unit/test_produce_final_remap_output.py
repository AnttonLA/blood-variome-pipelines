import os
import polars as pl
from produce_final_remap_output import produce_final_remap_output, filter_remap_output_file_by_biotype


def test_produce_final_remap_output():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    remap_lookup_outputs_dir = os.path.join(tests_dir, "test_data/dummy_remap_lookup_outputs/")
    output_file = os.path.join(tests_dir, "tmp/tmp_test_remap_final_output.tsv")
    produce_final_remap_output(remap_lookup_outputs_dir, output_file)

    # Ensure the file has been created
    assert os.path.isfile(output_file)

    # Ensure the file has correct columns (chr, pos, study_accession, transcription_factor, biotype, distance_to_peak)
    df = pl.read_csv(output_file, separator="\t", has_header=True)
    assert df.columns == ["chr", "pos", "study_accession", "transcription_factor", "biotype", "distance_to_peak"]

    # Compare "tmp/tmp_test_remap_final_output.tsv" and "test_data/dummy_remap_final_output.tsv"
    expected_df = pl.read_csv(os.path.join(tests_dir, "test_data/dummy_remap_final_output.tsv"), separator="\t",
                              has_header=True)
    assert df.frame_equal(expected_df)

    # Delete file once test is done
    os.remove(output_file)


def test_filter_remap_output_file_by_biotype():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    input_file = os.path.join(tests_dir, "test_data/dummy_remap_final_output.tsv")
    filter_file = os.path.join(tests_dir, "test_data/dummy_biotype_list.txt")
    reference_output_file = os.path.join(tests_dir, "test_data/dummy_remap_final_output_filtered.tsv")
    output_file = os.path.join(tests_dir, "tmp/tmp_test_remap_final_output_filtered.tsv")

    filter_remap_output_file_by_biotype(input_file, filter_file, output_file)

    # Check that the output and the reference are the same
    df = pl.read_csv(output_file, separator="\t", has_header=True)
    expected_df = pl.read_csv(reference_output_file, separator="\t", has_header=True)
    assert df.frame_equal(expected_df)

    # Delete file once test is done
    os.remove(output_file)
