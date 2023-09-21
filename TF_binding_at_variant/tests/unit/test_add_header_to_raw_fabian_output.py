from fabian.add_header_to_raw_fabian_output import add_header_to_raw_fabian_output


def test_add_header_to_raw_fabian_output():
    file = "test_data/dummy_raw_fabian_output.tsv"
    out_df = add_header_to_raw_fabian_output(file)

    assert out_df.columns == ['variant', 'tf', 'model_id', 'database', 'wt_score', 'mt_score', 'start_wt', 'end_wt',
                              'start_mt', 'end_mt', 'strand_wt', 'strand_mt', 'prediction', 'score']
