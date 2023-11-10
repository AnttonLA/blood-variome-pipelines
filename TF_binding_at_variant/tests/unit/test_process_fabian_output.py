import os
from process_fabian_output import process_fabian_output_data, add_header_to_raw_fabian_output_data,\
    process_fabian_output_table


def test_add_header_to_raw_fabian_output_data():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    file = os.path.join(tests_dir, "test_data/dummy_FABIAN_OUTPUT_data.tsv")
    out_df = add_header_to_raw_fabian_output_data(file)

    assert out_df.columns == ['variant', 'tf', 'model_id', 'database', 'wt_score', 'mt_score', 'start_wt', 'end_wt',
                              'start_mt', 'end_mt', 'strand_wt', 'strand_mt', 'prediction', 'score']


def test_process_fabian_output_data():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    file = os.path.join(tests_dir, "test_data/dummy_FABIAN_OUTPUT_data.tsv")
    map_file = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv.map")

    process_fabian_output_data(file, map_file, os.path.join(tests_dir, "tmp/dummy_FABIAN_OUTPUT_data.tsv.processed"))

    assert os.path.isfile(os.path.join(tests_dir, "tmp/dummy_FABIAN_OUTPUT_data.tsv.processed"))

    # Delete file once test is done
    os.remove(os.path.join(tests_dir, "tmp/dummy_FABIAN_OUTPUT_data.tsv.processed"))


def test_process_fabian_output_table():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    file = os.path.join(tests_dir, "test_data/dummy_FABIAN_OUTPUT_table.tsv")
    map_file = os.path.join(tests_dir, "test_data/dummy_variant_list.tsv.map")

    output_file = os.path.join(tests_dir, "tmp/test_process_fabian_output_table.tsv")
    process_fabian_output_table(file, map_file, output_file)

    assert os.path.isfile(output_file)

    # Delete file once test is done
    os.remove(output_file)

