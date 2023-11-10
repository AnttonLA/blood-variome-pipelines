import os
import polars as pl


def add_header_to_raw_fabian_output_data(headerless_file: str) -> pl.DataFrame:
    """
    Add the column names to the FABIAN-Variant output data file. This is necessary since the raw FABIAN output does not
    include the column names.
    Names taken from: https://www.genecascade.org/fabian/documentation#download-format
    Note that this function returns a DataFrame and does not write to file.

    :param headerless_file: Unmodified Fabian output file

    :return: Pandas DataFrame with added header
    """
    # Read the input TSV file without a header
    df = pl.read_csv(headerless_file, separator='\t', has_header=False)

    # Check if the number of columns matches the header columns
    num_columns = df.width
    expected_columns = 14

    if num_columns != expected_columns:
        raise ValueError(f"The input file should have {str(expected_columns)} columns, but has {num_columns} instead.")

    # Rename the columns
    rename_map = {"column_1": "variant", "column_2": "tf", "column_3": "model_id", "column_4": "database",
                  "column_5": "wt_score", "column_6": "mt_score", "column_7": "start_wt", "column_8": "end_wt",
                  "column_9": "start_mt", "column_10": "end_mt", "column_11": "strand_wt", "column_12": "strand_mt",
                  "column_13": "prediction", "column_14": "score"}
    df = df.rename(rename_map)

    return df


def process_fabian_output_data(fabian_output_file: str, map_file: str, output_file: str,
                               s_threshold: float = 0.2) -> None:
    """
    Process the raw FABIAN-Variant output file and produce a file with the following columns:
    - ID
    - Chrom
    - Pos
    - OA
    - EA
    - variant
    - tf
    - prediction
    - score


    :param fabian_output_file: Raw FABIAN-Variant output file
    :param map_file: Map file produced by variant_list_to_fabian_input_vcf()
    :param s_threshold: Threshold for the score column. Entries with an absolute score lower than this will be filtered
     out.
    :param output_file: Name of the output file that will be created

    :return:
    """
    # Check if the input file exists
    if not os.path.isfile(fabian_output_file) or os.path.getsize(fabian_output_file) == 0:
        raise ValueError(f"Input file {fabian_output_file} does not exist or is empty")

    # Add the header to the raw FABIAN-Variant output file
    df = add_header_to_raw_fabian_output_data(fabian_output_file)

    # Remove the character "*" in the "score" column
    df = df.with_columns(pl.col("score").str.replace("\*", ""))

    # Remove the suffix from the variant column
    df = df.with_columns(pl.col("variant").str.split(".").map_elements(lambda s: s[0]).alias("variant"))

    # Read the map file
    map_df = pl.read_csv(map_file, separator='\t', has_header=True)

    # Join the two DataFrames
    df = df.join(map_df, left_on="variant", right_on="Chrom:PosOA>EA", how="inner")

    # Rename "variant" to "Chrom:PosOA>EA"
    df = df.rename({"variant": "Chrom:PosOA>EA"})

    # Keep only the necessary columns
    df = df.select(["ID", "Chrom", "Pos", "OA", "EA", "Chrom:PosOA>EA", "tf", "prediction", "score"])

    # Filter out entries with a score lower than the threshold
    df = df.with_columns(pl.col("score").cast(pl.Float64))
    df = df.filter(pl.col("score").abs() >= s_threshold)

    # Write the output file
    df.write_csv(output_file, separator='\t', has_header=True)


def process_fabian_output_table(fabian_output_file: str, map_file: str, output_file: str,
                                s_threshold: float = 0.2) -> None:
    """
    Process the raw FABIAN-Variant output file and replace the column names with the proper variant IDs.
    This requires to melt the table and then join the rows by the FABIAN input format (Chrom:PosOA>EA) entries with the
    map file.
    The output file will have the following columns: ID, Chrom, Pos, OA, EA, Chrom:PosOA>EA, TF, score

    :param fabian_output_file:
    :param map_file:
    :param output_file:
    :param s_threshold:
    :return:
    """
    # Check if the input files exists and are non-empty
    if not os.path.isfile(fabian_output_file) or os.path.getsize(fabian_output_file) == 0:
        raise ValueError(f"Input file {fabian_output_file} does not exist or is empty")
    if not os.path.isfile(map_file) or os.path.getsize(map_file) == 0:
        raise ValueError(f"Input file {map_file} does not exist or is empty")

    df = pl.read_csv(fabian_output_file, separator='\t', has_header=True, infer_schema_length=0)
    # infer_schema_length=0 is necessary to avoid parsing issues.
    # See https://stackoverflow.com/questions/71106690/polars-specify-dtypes-for-all-columns-at-once-in-read-csv

    # Rename the first column to "TF"
    df = df.rename({"": "TF"})
    map_df = pl.read_csv(map_file, separator='\t', has_header=True)

    # Transpose the DataFrame and make the first row the header
    df = df.transpose(include_header=True, header_name="variant", column_names=df.get_column("TF").to_list())

    # Remove the first row, since it is the same as the header
    df = df[1:]

    # Remove the suffix from the variant column
    df = df.with_columns(pl.col("variant").str.split(".").map_elements(lambda s: s[0]).alias("variant"))

    # Melt the DataFrame. This is the most important step!
    out_df = df.melt(id_vars=["variant"], variable_name="TF", value_name="score")

    # Remove the character "*" in the "score" column
    out_df = out_df.with_columns(pl.col("score").str.replace("\*", ""))

    # Cast score to float and filter by score
    out_df = out_df.with_columns(pl.col("score").cast(pl.Float64))
    out_df = out_df.filter(pl.col("score").abs() >= s_threshold)

    # Join the two DataFrames
    out_df = out_df.join(map_df, left_on="variant", right_on="Chrom:PosOA>EA", how="inner")

    # Rename "variant" to "Chrom:PosOA>EA" and reorder the columns
    out_df = out_df.rename({"variant": "Chrom:PosOA>EA"})
    out_df = out_df.select(["ID", "Chrom", "Pos", "OA", "EA", "Chrom:PosOA>EA", "TF", "score"])

    # Remove the chr prefix from the Chrom column, cast it to int and sort by Chrom and Pos
    out_df = out_df.with_columns(pl.col("Chrom").str.replace("chr", "").cast(pl.Int64))
    out_df = out_df.sort(by=["Chrom", "Pos"])

    # Write to file
    out_df.write_csv(output_file, separator='\t', has_header=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process the raw FABIAN-Variant output file and produce a file with "
                                                 "the following columns: ID, Chrom, Pos, OA, EA, variant, tf, "
                                                 "prediction, score")
    parser.add_argument("-t", "--fabian_output_table_file", type=str, help="Raw FABIAN-Variant output TABLE file. This "
                                                                           "contains the summarized results.")
    parser.add_argument("-d", "--fabian_output_data_file", type=str, help="Raw FABIAN-Variant output DATA file. This "
                                                                          "contains the full results.")
    parser.add_argument("-m", "--map_file", required=True, help="Map file (ideally produced earlier by "
                                                                "variant_list_to_fabian_input_vcf()")
    parser.add_argument("-s", "--s_threshold", type=float, default=0.2, help="Threshold for the score column. Entries "
                                                                             "with an absolute score lower than this "
                                                                             "will be filtered out.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory where the output file will be created")

    args = parser.parse_args()

    if args.fabian_output_table_file is not None:
        output_file = os.path.join(args.output_dir, "fabian_output_table.processed")
        process_fabian_output_table(args.fabian_output_table_file, args.map_file, output_file, args.s_threshold)

    if args.fabian_output_data_file is not None:
        output_file = os.path.join(args.output_dir, "fabian_output_data.processed")
        process_fabian_output_data(args.fabian_output_data_file, args.map_file, output_file, args.s_threshold)
