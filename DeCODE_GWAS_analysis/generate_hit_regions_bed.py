import os
import sys
import polars as pl
import subprocess
import argparse

"""
This script is used to generate a .bed file that contains the regions where at least one GWAS hit is present.
Hits that are in within 1Mb of one other are merged into a single region, no matter the associated trait.
The output file also includes the names of the traits associated to the hits in each regions, and the smallest p-value
among the hits in the region.

Note that it uses 'bedtools'. You will need to give the path to wherever your bedtools is for the subprocess call in the
script to work.
"""

parser = argparse.ArgumentParser(description="Make a BED file listing genomic regions containing GWAS hits noted in the"
                                             " summary statistics files produced by 'extract_variants_by_pval.py'.")
parser.add_argument("folder_path",
                    metavar="FOLDER_PATH",
                    help="Path to the directory containing the summary stats files that contain the GWAS hits.")
parser.add_argument("-b", "--bedtools_path",  # TODO: IMPLEMENT THIS!
                    metavar="BEDTOOLS_PATH",
                    help="Path to the bedtools executable. If none is given, the script will assume that bedtools is "
                         "installed and available in the PATH.")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output file will be written to.")

args = parser.parse_args()


def print_status(percent):
    """Prints a status bar to the console. Used to show progress of the script when we iterate through files."""
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()


# Create empty output df, specify datatypes of each column.
bed_df = pl.DataFrame({"chrom": [],
                       "chromStart": [],
                       "chromEnd": [],
                       "leadSnp_pos": [],
                       "leadSnp_pval": [],
                       "phenotypes": []},
                      columns=[("chrom", pl.Int64), ("chromStart", pl.Int64), ("chromEnd", pl.Int64),
                               ("leadSnp_pos", pl.Int64), ("leadSnp_pval", pl.Float64), ("phenotypes", pl.Utf8)])

print("Welcome!\nReading all summary stats files in the provided folder...")
for i, file in enumerate(os.listdir(args.folder_path)):
    df = pl.read_csv(args.folder_path + file, sep="\t", columns=["pval", "chromosome", "position", "phenotype"])
    for row in df.rows():
        # row is a tuple(?) with this shape: ("pval", "chromosome", "position", "phenotype")
        bed_df = bed_df.extend(pl.DataFrame({"chrom": [int(row[1])],
                                             "chromStart": [max(2, row[2] - 10 ** 6)],
                                             "chromEnd": [row[2] + 10 ** 6],
                                             "leadSnp_pos": [row[2]],
                                             "leadSnp_pval": [row[0]],
                                             "phenotypes": [row[3]]}))

    if not (i % 10):  # Print status every 10 files
        percentage = 100 * (i / len(os.listdir(args.folder_path)))
        print_status(percentage)
    # if i >= 7:  # TODO: comment this out
    #     break

bed_df = bed_df.sort([pl.col("chrom"), pl.col("chromStart")])

# Add ID column to keep track of each row later.
id_list = [num for num in range(len(bed_df))]
bed_df.insert_at_idx(6, pl.Series(name="ID", values=id_list))

# Add column 'id_pos_pval_pheno' to keep track of rows with the same position and/or pval, which bedtools will merge.
bed_df = bed_df.with_column((pl.col("ID").cast(pl.Utf8) + ":" +
                             pl.col("chrom").cast(pl.Utf8) + ":" +
                             pl.col("leadSnp_pos").cast(pl.Utf8) + "&" +
                             pl.col("leadSnp_pval").cast(pl.Utf8) + "&" +
                             pl.col("phenotypes")).alias("id_pos_pval_pheno"))
# Generate bed file with a 1Mb window around each hit, up and down. The overlapping ranges will be merged with bedtools.
bed_df.write_csv(args.output_filepath + "/variant_regions_step1.bed", sep="\t", has_header=True)

########################################################################################################################
# Use bedtools to collapse the overlapping ranges while keeping the phenotype names
# See documentation: https://bedtools.readthedocs.io/en/latest/content/tools/merge.html?highlight=merge
print("Complete!\nMerging all overlapping regions with 'bedtools'...")

with open(args.output_filepath + "/variant_regions_step2.bed", "w") as f:
    subprocess.run(["/home/antton/Programs/bedtools2/bin/bedtools", "merge", "-i",
                    f"{args.output_filepath}/variant_regions_step1.bed", "-header", "-c", "4,5,6,7,8", "-o",
                    "distinct"],
                   check=True, stdout=f)

# Read the collapsed bed file and convert it to a polars dataframe
collapsed_bed_df = pl.read_csv(args.output_filepath + "/variant_regions_step2.bed", sep="\t",
                               columns=["chrom", "chromStart", "chromEnd", "leadSnp_pos", "leadSnp_pval",
                                        "phenotypes", "ID", "id_pos_pval_pheno"])

# Extract position and p-value lists from the 'id_pos_pval_pheno' column
collapsed_bed_df = collapsed_bed_df.with_column((pl.col("id_pos_pval_pheno").apply(
    lambda x: ",".join([y.split("&")[0].split(":")[-1] for y in x.split(",")])).alias("pos_list")))
collapsed_bed_df = collapsed_bed_df.with_column(
    (pl.col("id_pos_pval_pheno").apply(lambda x: ",".join([y.split("&")[1] for y in x.split(",")])).alias("pval_list")))

# Add two new columns: total number of hits in a region, and number of phenotypes linked to those hits.
collapsed_bed_df.insert_at_idx(5, collapsed_bed_df.select("ID")
                               .apply(lambda x: len(str(x).split(',')) - 1).to_series().alias("num_included_variants"))
collapsed_bed_df.insert_at_idx(6, collapsed_bed_df.select("phenotypes")
                               .apply(lambda x: len(str(x).split(',')) - 1).to_series().alias("num_phenotypes"))


########################################################################################################################
# Keep only the lead SNPs position and pval info for each region, instead of the full lists given by bedtools.
def find_row_min_index(comma_separated_str) -> int:
    """Returns the index of the minimum value in a string of comma-separated numbers."""
    str_as_str_list = comma_separated_str.split(',')
    str_as_float_list = [float(x) for x in str_as_str_list]
    min_value = min(str_as_float_list)
    return str_as_float_list.index(min_value)


def num_elements_in_str(comma_separated_str) -> int:
    """Returns the number of elements in a string of comma-separated numbers."""
    return len(comma_separated_str.split(','))


# Create subset of the dataframe with only the columns we need to find the lead SNPs
leadSnp_df = collapsed_bed_df.select(["ID", "pos_list", "pval_list"])

# Check concordant number of elements in each list:
leadSnp_df = leadSnp_df.with_column(
    (pl.col("ID").apply(num_elements_in_str) == pl.col("pval_list").apply(num_elements_in_str)).alias("check"))
if not leadSnp_df.select(pl.col("check").all())[0, 0]:  # True if all elements are True
    raise ValueError("Number of elements in columns 'ID' and 'pval_list' columns do not match!")

# Find the index of the lowest p-value in each row
leadSnp_df = leadSnp_df.with_column(pl.col("pval_list").apply(find_row_min_index).alias("leadSnp_index"))

# Extract the lead SNP position and p-value from the lists
# Create columns and convert contents to list so that arr.get() can be used
leadSnp_df = leadSnp_df.with_columns(
    [pl.col("pval_list").apply(lambda x: x.split(',')).alias("leadSnp_pval"),
     pl.col("pos_list").apply(lambda x: x.split(',')).alias("leadSnp_pos")])
# Extract values at desired index
leadSnp_df = leadSnp_df.with_columns(
    [pl.col("leadSnp_pval").arr.get(pl.col("leadSnp_index")).alias("leadSnp_pval"),
     pl.col("leadSnp_pos").arr.get(pl.col("leadSnp_index")).alias("leadSnp_pos")])

# Replace the original columns with the new ones
collapsed_bed_df.replace("leadSnp_pos", leadSnp_df.select("leadSnp_pos").to_series())
collapsed_bed_df.replace("leadSnp_pval", leadSnp_df.select("leadSnp_pval").to_series())

# Remove the columns we don't need anymore
collapsed_bed_df = collapsed_bed_df.drop(["ID", "id_pos_pval_pheno", "pos_list", "pval_list"])


# print("collapsed_bed_df: ", collapsed_bed_df)
# print(tuple(zip(collapsed_bed_df.columns, collapsed_bed_df.row(2))))
# sys.exit()


output_file_name = "variant_regions.bed"
collapsed_bed_df.write_csv(args.output_filepath + '/' + output_file_name, sep="\t", has_header=True)
print(f"Done! Final regions saved to '{args.output_filepath + '/' + output_file_name}'")

# Remove the intermediate files
os.remove(args.output_filepath + "/variant_regions_step1.bed")
os.remove(args.output_filepath + "/variant_regions_step2.bed")
