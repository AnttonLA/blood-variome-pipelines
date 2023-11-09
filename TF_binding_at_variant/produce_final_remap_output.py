import os
import sys
import polars as pl
import argparse

"""
This file will combine every individual file in the 'remap_lookup_outputs/' directory into a single file.
The final file will have the following columns:
chr | pos | study_accession | transcription_factor | biotype | distance_to_snp

The script will also offer the chance to filter the entries by biotype.

This script is expected to be run after 'remap_lookup_for_full_snplist.py'.
"""


def produce_final_remap_output(interim_file_dir: str, output_file: str) -> None:
    """
    This function will combine every individual file in the 'remap_lookup_outputs/' directory into a single file.
    The files have the following naming convention: remap_studies_<chr:pos>.txt
    And the following columns:

    study_accession | transcription_factor | biotype | distance_to_peak

    The output file will combine all of them into a single file and have the following columns:

    chr | pos | study_accession | transcription_factor | biotype | distance_to_snp

    :param interim_file_dir: Directory where the individual files are stored. Usually named 'remap_lookup_outputs/'.
    :param output_file: Name of the output file.
    :return:
    """
    # Ensure interim_file_dir exists and is not empty
    if not os.path.isdir(interim_file_dir):
        raise ValueError("Directory does not exist")
    if len(os.listdir(interim_file_dir)) == 0:
        raise ValueError("Directory is empty")

    files_list = os.listdir(interim_file_dir)  # Get the name of every file in the remap_lookup_outputs/ directory
    # Empty DataFrame with the correct schema
    df = pl.DataFrame(schema=[("chr", pl.Int32), ("pos", pl.Int32), ("study_accession", pl.Utf8),
                              ("transcription_factor", pl.Utf8), ("biotype", pl.Utf8), ("distance_to_peak", pl.Int32)])
    # Fill up the DataFrame with the contents of each file
    for file in files_list:
        # Get chr:pos from the filename
        chr_pos = file.replace("remap_studies_", "").replace(".txt", "")
        chr_int = int(chr_pos.split(":")[0])
        pos_int = int(chr_pos.split(":")[1])

        tmp_df = pl.read_csv(os.path.join(interim_file_dir, file), separator="\t", has_header=True,
                             dtypes={"study_accession": pl.Utf8, "transcription_factor": pl.Utf8, "biotype": pl.Utf8,
                                     "distance_to_peak": pl.Int32})
        tmp_df = tmp_df.with_columns([pl.lit(chr_int).alias("chr"), pl.lit(pos_int).alias("pos")])
        # Reorder columns
        tmp_df = tmp_df.select(["chr", "pos", "study_accession", "transcription_factor", "biotype", "distance_to_peak"])
        df = df.vstack(tmp_df)

    # Sort by chr and pos
    df = df.sort(["chr", "pos"])

    # Write to file
    df.write_csv(output_file, separator="\t")


def filter_remap_output_file_by_biotype(full_file: str, biotypes: list[str], output_file: str) -> None:
    """
    This function will filter the final output file by biotype.

    :param full_file: Path to the full output file
    :param biotypes: List of biotypes to filter by
    :param output_file: Name of the output file
    :return:
    """
    # Ensure full_file exists and is not empty
    if not os.path.isfile(full_file):
        raise ValueError("File does not exist")
    if os.path.getsize(full_file) == 0:
        raise ValueError("File is empty")

    # If biotypes is empty, the output file will be the same as the input file. Issue a warning and return.
    if len(biotypes) == 0:
        sys.stderr.write(f"\nWARNING: No biotype filters were provided. Output file {output_file} will be identical"
                         f"to full ReMap lookup file {full_file}.\n")
    else:
        # Read the full file into a DataFrame
        df = pl.read_csv(full_file, separator="\t", has_header=True,
                         dtypes={"chr": pl.Int32, "pos": pl.Int32, "study_accession": pl.Utf8,
                                 "transcription_factor": pl.Utf8, "biotype": pl.Utf8, "distance_to_peak": pl.Int32})
        # Filter the DataFrame by biotype
        df = df.filter(pl.col("biotype").is_in(biotypes))
        # Write to file
        df.write_csv(output_file, separator="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Look up the ReMap database and produce an output file with all the TFs that bind"
                                     "over the variants on interest. It is possible to filter the output by biotype.")
    parser.add_argument("-d", "--interim_file_dir", type=str, help="Directory where the individual files are stored. "
                                                                   "Usually named 'remap_lookup_outputs/'.")
    parser.add_argument("-f", "--filter", nargs="+", help="OPTIONAL. List of biotypes to filter by. If none is "
                                                          "specified, the output file will be the "
                                                          "complete output of the lookup.")
    parser.add_argument("-o", "--output_file", required=True, help="Name of the output file.")

    args = parser.parse_args()

    # Check that the interim file directory exists and is not empty
    if not os.path.isdir(args.interim_file_dir):
        raise ValueError("Directory does not exist")
    if len(os.listdir(args.interim_file_dir)) == 0:
        raise ValueError("Directory is empty")

    produce_final_remap_output(args.interim_file_dir, args.output_file)
    sys.stdout.write(f"\nWrote full ReMap lookup results to: {args.output_file}\n")

    if args.filter:

        # Check that the filter is a list
        if not isinstance(args.filter, list):
            raise ValueError("Filter must be a list")

        sys.stdout.write(f"\nFiltered file requested. Creating second file with only the biotypes: {args.filter}\n")

        filtered_file_name = args.output_file.replace(".tsv", "_filtered.tsv")
        filter_remap_output_file_by_biotype(args.output_file, args.filter, filtered_file_name)
        sys.stdout.write(f"\nWrote filtered ReMap lookup results to: {filtered_file_name}\n")
