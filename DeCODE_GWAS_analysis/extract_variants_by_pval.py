import os
import sys

import polars as pl
from scipy.stats import chi2
import argparse

"""
This script takes a 'variant_info.txt' file and a number of GWAS output .rse files (ID, beta and chi-square columns).
It then filters the output tables by a user specified p-value threshold.
The resulting entries, together with additional info about the variants, are written to a new file.
"""

parser = argparse.ArgumentParser(description="Extract the entries with p-values lower than a specified threshold from "
                                             "a set of GWAS output files.")
parser.add_argument("variant_info_file",
                    metavar="FILEPATH",
                    help="Path to the variant info file.")
parser.add_argument("files_filepath",
                    metavar="FOLDER_PATH",
                    help="Path to the directory containing the GWAS output files.")
parser.add_argument("-p", "--pval_thresh", default=6,
                    help="Negative logarithm of the p-value. This will be the exponent of the desired p-value "
                         "threshold. E.g. 6 for p=1e-6. Only variants with p-values smaller than the corresponding "
                         "threshold will be included in the output tables. Default: 6")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output files will be written to.")

args = parser.parse_args()

# TODO: Input validation


def print_status(percent):
    """Prints a status bar to the console. Used to show progress of the script when we iterate through files."""
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()


# Read variant info file, then merge it with the GWAS output files
var_info_df = pl.read_csv(args.variant_info_file, sep='\t', columns=["ID", "Marker", "OA", "EA", "EAF", "Info"])

pval_threshold = 10 ** -float(args.pval_thresh)
chi2_threshold = chi2.isf(pval_threshold, 1)


for i, gwas_file in enumerate(os.listdir(args.files_filepath)):

    # Progress bar
    percentage = int(i / len(os.listdir(args.files_filepath)) * 100)
    print_status(percentage)

    if gwas_file.endswith('.res'):  # TODO: Maybe allow more extensions?
        gwas_df = pl.read_csv(args.files_filepath + '/' + gwas_file, has_header=False,
                              new_columns=["ID", "beta", "chi2"], sep=" ")
        gwas_df = gwas_df.filter(gwas_df["chi2"] > chi2_threshold)
        if gwas_df.is_empty():  # If there are no hits, write an empty file to appease Snakemake # TODO: sure it works?
            mock_df = pl.DataFrame(columns=["ID", "beta", "chi2", "pval", "Marker", "chromosome", "position", "OA",
                                            "EA", "EAF", "Info", "phenotype"])
            phenotype = '_'.join(
                gwas_file.rstrip('.txt').split('_')[4:-3])  # Phenotype name is only given by the file name
            mock_df.write_csv(args.output_filepath + '/' + phenotype + '.txt', sep='\t', has_header=True)
            continue

        # Add 'pval' column
        gwas_df.insert_at_idx(3, pl.Series("pval", [chi2.sf(x, 1) for x in gwas_df["chi2"]]))

        # Merge with variant info
        gwas_df = gwas_df.join(var_info_df, on="ID", how="left")

        # Add chromosome and position columns, and sort the table with them.
        # We replace 'chrX' with 'chr23' to allow sorting. We assume marker format to be chr<chr>:<pos>
        gwas_df = gwas_df.with_column(pl.col("Marker").apply(lambda x: "chr23"+x[4:] if x.startswith("chrX") else x))
        gwas_df.insert_at_idx(5, pl.Series("chromosome", [int(x.split(':')[0][3:]) for x in gwas_df["Marker"]]))
        gwas_df.insert_at_idx(6, pl.Series("position", [int(x.split(':')[1]) for x in gwas_df["Marker"]]))
        gwas_df = gwas_df.sort(["chromosome", "position"])

        # Add 'phenotype' column
        phenotype = '_'.join(gwas_file.rstrip('.txt').split('_')[4:-3])  # Phenotype name is only given by the file name
        gwas_df.hstack([pl.Series("phenotype", [phenotype] * gwas_df.height)], in_place=True)

        # Write to file
        gwas_df.write_csv(args.output_filepath + '/' + phenotype + '.txt', sep='\t', has_header=True)
