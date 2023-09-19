from extract_remapdb_studies import extract_studies_for_single_snp


def test_extract_studies_for_single_snp():
    extract_studies_for_single_snp("chr1:23935190", "../test_data/dummy_remap_file.bed.gz", "../test_data/tmp_dir",
                                   "../test_data/output.tsv", True)
    assert False
