import polars as pl
import os
import argparse
import sys


def extract_studies_for_single_snp(chr_pos: str, remap_file: str, tmp_dir: str, output: str, verbose: bool = False) -> None:
    """
    This function will take the position of a SNP ('chr_pos') and produce a tabix query to the ReMap metadata file
    (filepath stored in 'remap_file'). It will then process the output of the tabix query to produce a table with the
    following columns: study_accession | transcription_factor | biotype | distance_to_snp. The table will be written to
    a text file (filepath stored in 'output').

    :param chr_pos: Position of the SNP in the format <chr>:<pos>
    :param remap_file: Path to ReMap BED file. Requires tabix index file.
    :param tmp_dir: Directory to store temporary files (i.e. tabix "slices" of the remap metadata file for the queried
     position). If none is specified, temporary files will be stored in the current working directory.
    :param output: Output file.
    :param verbose: If True, print progress to stdout.
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

    if verbose:
        sys.stdout.write(f"Requested position: <{snp_full_pos}>\n")

    # Check if ReMap file exists, it is bgzipped, and it has a tabix index file
    if not os.path.isfile(remap_file):
        raise ValueError("ReMap file does not exist")
    if not remap_file.endswith(".gz"):
        raise ValueError("ReMap file must be bgzipped")
    if not os.path.isfile(remap_file + ".tbi"):
        raise ValueError("No tabix index file could be found for the ReMap file.")

    # Start processing the tabix query output
    # Compose tabix query
    command = f"tabix {remap_file} {snp_full_pos}"

    # Run command and store output. Requires tabix.
    query_output = os.popen(command).read()

    # Check if the output is empty. If so, give a warning, write an empty file, and return None to quit early.
    if query_output == "":
        sys.stderr.write(f"WARNING in extract_studies_for_single_snp : "
                         f"No ReMap entries found for position {snp_full_pos}\n")
        if verbose:
            sys.stderr.write(f"Wrote empty file: {output}\n")
        # Write empty file
        with open(output, "w") as f:
            f.write("study_accession\ttranscription_factor\tbiotype\tdistance_to_peak\n")
        return None

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

    df = pl.read_csv(tmp_file_path, sep="\t", has_header=False,
                     new_columns=["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
                                  "itemRgb"],
                     dtype={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "name": pl.Utf8, "score": pl.Float64,
                            "strand": pl.Utf8, "thickStart": pl.Int64, "thickEnd": pl.Int64, "itemRgb": pl.Utf8})

    # Distance of SNP to peak. We'll use this to sort entries later
    # It'd probably be fine to calculate distance to thickStart, but we'll calculate the center of the peak to be safe
    df = df.with_column(pl.struct(['thickStart', 'thickEnd'])
                        .apply(lambda s: int((s['thickStart'] + s['thickEnd']) / 2)).cast(pl.Int64)
                        .alias("thickCenter"))

    # Calculate distance to SNP pos (stored in variable 'pos')
    df = df.with_column(pl.col("thickCenter")
                        .apply(lambda s: abs(int(s) - int(pos))).cast(pl.Int64)
                        .alias("distance_to_peak"))

    # Split the name column by '.' and make new columns out of the three entries
    df = df.with_columns([pl.col("name").str.split(".").apply(lambda s: s[0]).alias("study_accession"),
                          pl.col("name").str.split(".").apply(lambda s: s[1]).alias("transcription_factor"),
                          pl.col("name").str.split(".").apply(lambda s: s[2]).alias("biotype")])

    out_df = df.select(["study_accession", "transcription_factor", "biotype", "distance_to_peak"])

    # Sort by distance from SNP to ChIP-seq peak.
    out_df = out_df.sort("distance_to_peak")

    # Write to file
    out_df.write_csv(output, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chr_pos", help="Position of the SNP in the format <chr>:<pos>")
    parser.add_argument("-r", "--remap_file", required=True,
                        help="Path to ReMap BED file. Requires tabix index file.")
    parser.add_argument("-t", "--tmp_dir", help="OPTIONAL. Directory to store temporary files. If none is specified, "
                                                "temporary files will be stored in the current working directory.")
    parser.add_argument("-o", "--output", required=True, help="Output file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output.")

    args = parser.parse_args()

    extract_studies_for_single_snp(args.chr_pos, args.remap_file, args.tmp_dir, args.output, args.verbose)
