import os

import polars as pl
from scipy.stats import chi2
import argparse

"""
This script takes a 'variant_info.txt' file and a GWAS output .rse file (ID, beta and chi-square columns), and combines
them to produce a full summary statistics file. The columns of the output file are: ID, beta, chi2, pval, Marker,
chromosome, position, OA, EA, EAF, Info and	phenotype.
Note that this script might take a while to run, as it needs to do calculations for each row of the GWAS output file.
"""

parser = argparse.ArgumentParser(description="Extract the contents of a 'variant_info.txt' and a '.res' file into a"
                                             " full summary statistics file for a single trait.")
parser.add_argument("variant_info_file",
                    metavar="FILEPATH",
                    help="Path to the variant info file.")
parser.add_argument("files_filepath",
                    metavar="FOLDER_PATH",
                    help="Path to the directory containing the GWAS output files.")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output file will be written to.")

args = parser.parse_args()

# TODO: Input validation

# Read variant info file, then merge it with the GWAS output file
var_info_df = pl.read_csv(args.variant_info_file, sep='\t', columns=["ID", "Marker", "OA", "EA", "EAF", "Info"])

# TODO: Add the option to select the phenotype we want
gwas_file_name = os.listdir(args.files_filepath)[0]  # We take the first file. It doesn't matter which one we take.
print("Commencing generation of full summary statistics file. This can take several minutes.\n"
      f"The selected trait was {gwas_file_name}")

gwas_df = pl.read_csv(args.files_filepath + '/' + gwas_file_name, has_header=False,
                      new_columns=["ID", "beta", "chi2"], sep=" ")

print("Finished loading GWAS output file. Adding p-values column...")
# Add 'pval' column
gwas_df.insert_at_idx(3, pl.Series("pval", [chi2.sf(x, 1) for x in gwas_df["chi2"]]))

print("Commencing merge with variant info file...")
# Merge with variant info
gwas_df = gwas_df.join(var_info_df, on="ID", how="left")

# Add chromosome and position columns, and sort the table with them.
# We replace 'chrX' with 'chr23' to allow sorting. We assume marker format to be chr<chr>:<pos>
gwas_df = gwas_df.with_column(pl.col("Marker").apply(lambda x: "chr23" + x[4:] if x.startswith("chrX") else x))
gwas_df.insert_at_idx(5, pl.Series("chromosome", [int(x.split(':')[0][3:]) for x in gwas_df["Marker"]]))
gwas_df.insert_at_idx(6, pl.Series("position", [int(x.split(':')[1]) for x in gwas_df["Marker"]]))
gwas_df = gwas_df.sort(["chromosome", "position"])

print("Almost done! Adding phenotype column...")
# Add 'phenotype' column
phenotype = '_'.join(gwas_file_name.rstrip('.txt').split('_')[4:-3])  # Phenotype name is only given by the file name
gwas_df.hstack([pl.Series("phenotype", [phenotype] * gwas_df.height)], in_place=True)

gwas_df.write_csv(args.output_filepath + '/template_manhattan.txt', sep='\t', has_header=True)
print("Done! Full summary statistics file written to " + args.output_filepath + '/template_manhattan.txt')
