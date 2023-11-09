import os
from extract_remapdb_studies import extract_studies_for_single_snp


def test_extract_studies_for_single_snp():
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.split(dir_of_file)[0]

    remap_file = os.path.join(tests_dir, "test_data/dummy_remap_file.bed.gz")
    tmp_dir = os.path.join(tests_dir, "tmp/")
    out_file = os.path.join(tests_dir, "test_data/output.tsv")

    extract_studies_for_single_snp("1:23939135", remap_file, tmp_dir, out_file, True)

    # Assert that the output file exists
    assert os.path.isfile(out_file)

    # Remove both the 'tabix_slices' directory, its contents, and the output file
    os.remove(out_file)
    for file in os.listdir(tmp_dir + "/tabix_slices/"):
        os.remove(tmp_dir + "/tabix_slices/" + file)
    os.rmdir(tmp_dir + "/tabix_slices")
