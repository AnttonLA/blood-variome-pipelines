import os
import sys
import polars as pl
import argparse

"""
This script combines all the summary statistics files produced by 'extract_variants_by_pval.py' into a single file.
"""

parser = argparse.ArgumentParser(description="Combine all summary statistics files produced by "
                                             "'extract_variants_by_pval.py' into one.")
parser.add_argument("files_filepath",
                    metavar="FOLDER_PATH",
                    help="Path to the directory containing the files to be merged.")
parser.add_argument("gwas_run_name",
                    metavar="GWAS_RUN_NAME",
                    help="Name of the GWAS run. This is used to name the output file.")
parser.add_argument("-p", "--pval_threshold", default="", type=str,
                    help="P-value threshold that was used when filtering variants. This is used to name the output"
                         " file. If non is given none will be added to the file name.")
parser.add_argument("-r", "--include_repeats", action="store_true",
                    help="by default, only the most significant phenotype will be included in the final file. If this "
                         "flag is present, all the significant associations for each variant will be included instead.")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output file will be written to.")

args = parser.parse_args()

full_hits_table_df = pl.DataFrame()

if not os.path.exists(args.files_filepath):
    raise ValueError("The given directory is empty.")

# Iterate through the files, concat contents into a df
num_all_tables = len(os.listdir(args.files_filepath))
for i, hit_table_file in enumerate(os.listdir(args.files_filepath)):
    if hit_table_file.endswith('.tsv') or hit_table_file.endswith('.txt'):
        # Table made from current file
        table_df = pl.read_csv(args.files_filepath + '/' + hit_table_file, sep="\t",
                               columns=["ID", "beta", "chi2", "pval", "Marker", "chromosome", "position", "OA", "EA",
                                        "EAF", "Info", "phenotype"])
        if table_df.is_empty():
            continue
        full_hits_table_df = pl.concat([full_hits_table_df, table_df])

        if (num_all_tables > 100) & (i % 100 == 0):  # Print every 100 if # of files > 100 files, every 10 if < 100
            print(f"Table  {i} / {num_all_tables}")
        elif (num_all_tables < 100) & (i % 10 == 0):
            print(f"Table  {i} / {num_all_tables}")
else:
    print(f"Table  {i+1} / {num_all_tables}\nDONE!")

if args.pval_threshold != "":
    pval_threshold_str = "_10E" + args.pval_threshold
else:
    pval_threshold_str = ""

if not args.include_repeats:  # Each variant only once
    # For each row with the same ID, keep the one with the lowest p-value. See README for more info.
    # This is the default behaviour.
    # TODO: ensure this works correctly
    full_hits_table_df = full_hits_table_df.sort(['pval']).groupby('ID').head(1)  # First is smallest after sorting
    full_hits_table_df = full_hits_table_df.sort(['chromosome', 'position'])

    hit_table_file_name = f'{args.gwas_run_name}_hits_only{pval_threshold_str}_one_pheno_only_per_variant.txt'
    print("Rows of the final table that contains all hits (counting each variant only once): ",
          len(full_hits_table_df["ID"]))

else:  # each variant can appear multiple times if it has multiple associations with different phenotypes
    full_hits_table_df = full_hits_table_df.sort(['chromosome', 'position'])
    hit_table_file_name = f'{args.gwas_run_name}_hits_only{pval_threshold_str}_all_phenotypes_per_variant.txt'
    print("Rows of the final table that contains all hits (same variant can be counted several times): ",
          len(full_hits_table_df["ID"]))

# Save the final table to file
full_hits_table_df.write_csv(args.output_filepath + '/' + hit_table_file_name, sep="\t")
