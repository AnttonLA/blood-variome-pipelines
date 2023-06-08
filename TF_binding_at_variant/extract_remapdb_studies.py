import polars as pl
import os
import argparse
import sys


def extract_studies_for_single_snp(chr_pos: str, remap_file: str, tmp_dir: str, output: str) -> None:
    """
    This function will take the position of a SNP ('chr_pos') and produce a tabix query to the ReMap metadata file
    (filepath stored in 'remap_file'). It will then process the output of the tabix query to produce a table with the
    following columns: study_accession | transcription_factor | biotype | distance_to_snp. The table will be written to
    a CSV file (filepath stored in 'output').

    :param chr_pos: Position of the SNP in the format <chr>:<pos>
    :param remap_file: Path to ReMap BED file. Requires tabix index file.
    :param tmp_dir: Directory to store temporary files. If none is specified, temporary files will be stored in the
    current working directory.
    :param output: Output file. Must be a CSV file.
    :return: None
    """

    # Input sanitation
    if ':' not in chr_pos:
        raise ValueError("SNP position must be in the format '<chr>:<position>'")

    if not chr_pos.split(":")[1].isdigit():
        raise ValueError("SNP position must be in the format '<chr>:<position>'")

    # Create 'snp_full_pos' variable
    chrom = chr_pos.split(":")[0]  # Chromosome of the SNP
    pos = chr_pos.split(":")[1]  # Position of the SNP
    snp_full_pos = f"chr{chrom}:{pos}-{pos}"  # Full position of SNP in 'chr<chrom>:<pos>-<pos>' format, for tabix query

    sys.stdout.write(f"Requested position: <{snp_full_pos}>\n")

    # Check if ReMap file exists, it is bgzipped, and it has a tabix index file
    if not os.path.isfile(remap_file):
        raise ValueError("ReMap file does not exist")
    if not remap_file.endswith(".gz"):
        raise ValueError("ReMap file must be bgzipped")
    if not os.path.isfile(remap_file + ".tbi"):
        raise ValueError("No tabix index file could be found for the ReMap file.")

    # Check if output file is a CSV file. If not, add extension '.csv' to output file name.
    if not output.endswith(".csv"):
        sys.stderr.write("WARNING: Output file must be a CSV file. Adding extension '.csv' to output file name.\n")
        output = output + ".csv"

    # Start processing the tabix query output
    # Compose tabix query
    command = f"tabix {remap_file} {snp_full_pos}"

    # Run command and store output. Requires tabix.
    query_output = os.popen(command).read()

    # Check if the output is empty
    if query_output == "":
        raise ValueError("No ReMap entries found for this SNP.")

    # Save the output of the tabix query to a file
    if not tmp_dir:  # If no tmp_dir is specified, use the current working directory
        tmp_dir = os.getcwd()
    # Check if a directory named 'tabix_slices' exists in the tmp_dir. If not, create it.
    if not os.path.isdir(f"{tmp_dir}/tabix_slices"):
        os.mkdir(f"{tmp_dir}/tabix_slices")
    tmp_file = f"remap2022_all_macs2_hg38_v1_0.{snp_full_pos}.bed"
    tmp_file_path = f"{tmp_dir}/tabix_slices/{tmp_file}"
    with open(tmp_file_path, "w") as f:
        f.write(query_output)

    df = pl.read_csv(tmp_file_path, sep="\t",
                     new_columns=["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
                                  "itemRgb"])

    # Distance of SNP to peak. We'll use this to sort entries later
    # It'd probably be fine to calculate distance to thickStart, but we'll calculate the center of the peak to be safe
    df = df.with_column(pl.struct(['thickStart', 'thickEnd'])
                        .apply(lambda s: int((s['thickStart'] + s['thickEnd']) / 2))
                        .alias("thickCenter"))

    # Calculate distance to SNP pos (stored in variable 'pos')
    df = df.with_column(pl.col("thickCenter").apply(lambda s: abs(int(s) - int(pos))).alias("distance_to_peak"))

    # Split the name column by '.' and make new columns out of the three entries
    df = df.with_columns([pl.col("name").str.split(".").apply(lambda s: s[0]).alias("study_accession"),
                          pl.col("name").str.split(".").apply(lambda s: s[1]).alias("transcription_factor"),
                          pl.col("name").str.split(".").apply(lambda s: s[2]).alias("biotype")])

    out_df = df.select(["study_accession", "transcription_factor", "biotype", "distance_to_peak"])

    # Sort by distance from SNP to ChIP-seq peak.
    out_df = out_df.sort("distance_to_peak")

    # Write to file
    out_df.write_csv(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chr_pos", help="Position of the SNP in the format <chr>:<pos>")
    parser.add_argument("-r", "--remap_file", required=True,
                        help="Path to ReMap BED file. Requires tabix index file.")
    parser.add_argument("-t", "--tmp_dir", help="OPTIONAL. Directory to store temporary files. If none is specified, "
                                                "temporary files will be stored in the current working directory.")
    parser.add_argument("-o", "--output", required=True, help="Output file. Must be a CSV file.")
    args = parser.parse_args()

    extract_studies_for_single_snp(args.chr_pos, args.remap_file, args.tmp_dir, args.output)
