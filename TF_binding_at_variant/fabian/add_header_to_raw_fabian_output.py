import polars as pl


def add_header_to_raw_fabian_output(input_file, output_file):
    """
    The default Fabian output does not inlcude the column names. This function will add the column names to the output.
    Names taken from: https://www.genecascade.org/fabian/documentation#download-format
    :param input_file:
    :param output_file:
    :return:
    """
    # Read the input TSV file without a header
    df = pl.read_csv(input_file, sep='\t', has_header=False)

    # Check if the number of columns matches the header columns
    num_columns = df.width
    expected_columns = 14

    if num_columns != expected_columns:
        print(f"Error: The input file should have {str(expected_columns)} columns.")
        return

    # Rename the columns
    rename_map = {"column_1": "variant", "column_2": "tf", "column_3": "model_id", "column_4": "database",
                  "column_5": "wt_score", "column_6": "mt_score", "column_7": "start_wt", "column_8": "end_wt",
                  "column_9": "start_mt", "column_10": "end_mt", "column_11": "strand_wt", "column_12": "strand_mt",
                "column_13": "prediction", "column_14": "score"}
    df = df.rename(rename_map)

    # Save the DataFrame to a new TSV file with the header
    df.write_csv(output_file, sep='\t', has_header=True)


# Example usage:
input_file = "/home/antton/Projects/CordBlood_GWAS/transcription_factor_lookup/fabian_output.tsv"
output_file = "/home/antton/Projects/CordBlood_GWAS/transcription_factor_lookup/fabian_output_h.tsv"
add_header_to_raw_fabian_output(input_file, output_file)
