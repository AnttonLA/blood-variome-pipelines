import os
import sys
import polars as pl
import subprocess
import argparse

"""
This script is used to generate a .bed-like file that contains the regions where at least one GWAS hit is present.
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
parser.add_argument("-s", "--interval_size",
                    metavar="INTERVAL_SIZE",
                    type=int,
                    default=1000000,
                    help="Size of the intervals in which the hits will be grouped. Default is 1000000 (1Mb) up and down"
                            " stream from each variant.")
parser.add_argument("-b", "--bedtools_path",
                    metavar="BEDTOOLS_PATH",
                    help="Path to the bedtools executable. If none is given, the script will assume that bedtools is "
                         "installed and available in the PATH.")
parser.add_argument("-t", "--tabix_path",
                    metavar="TABIX_PATH",
                    help="Path to the tabix executable. If none is given, the script will assume that tabix is "
                         "installed and available in the PATH.")
parser.add_argument("-g", "--gtf_file",
                    metavar="GTF_FILE",
                    help="Path to the GTF file that will be used to find the closest gene to each hit. If none is given"
                         ", the 'leadSnp_closest_gene' entry in the output file will be set to NA for every region.")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output file will be written to.")

args = parser.parse_args()

# Deal with the input arguments
if args.bedtools_path is None:
    bedtools_path = "bedtools"
if args.tabix_path is None:
    tabix_path = "tabix"


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
                                             "chromStart": [max(2, row[2] - args.interval_size)],
                                             "chromEnd": [row[2] + args.interval_size],
                                             "leadSnp_pos": [row[2]],
                                             "leadSnp_pval": [row[0]],
                                             "phenotypes": [row[3]]}))

    if not (i % 10):  # Print status every 10 files
        percentage = 100 * (i / len(os.listdir(args.folder_path)))
        print_status(percentage)

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
    subprocess.run([args.bedtools_path, "merge", "-i",
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


########################################################################################################################
# Include the closest gene to the lead SNP in each region
def find_closest_gene(snp_pos: str, tabix_output_list: list) -> str:
    """
    From a list of tabix output lines from the GENCODE v42 gtf file, returns the row pertaining to the closest gene to
    the given SNP position. This function is called by the 'get_closest_gene_tabix_output_row' function.

    :param snp_pos: string with the position of the SNP (just 'pos', no 'chr')
    :param tabix_output_list: list of strings with the output of a tabix query on the GENCODE v42 gtf file
    :return: string with the row referring to the closest gene to the given SNP position
    """
    closest_gene_entry_index = 0
    for index, line in enumerate(tabix_output_list):
        split_line = line.split("|")
        gene_start = int(split_line[3])
        distance = abs(gene_start - int(snp_pos))
        winning_distance = abs(int(tabix_output_list[closest_gene_entry_index].split("|")[3]) - int(snp_pos))
        if distance < winning_distance:
            closest_gene_entry_index = index

    return tabix_output_list[closest_gene_entry_index]


def get_closest_gene_tabix_output_row(x) -> str:
    """
    For each row, call tabix for that position, and return the first row of the output that corresponds to "gene".
    If there are none, return an empty string. Same if the tabix query has no output.
    This function is used with the 'apply' method of a polars dataframe
    """
    cromosome_num = x["chrom"]
    lead_snp_pos = x["leadSnp_pos"]
    tabix_query_range_start = str(int(lead_snp_pos) - 1000000)
    tabix_query_range_end = str(int(lead_snp_pos) + 1000000)
    out_str = subprocess.run([args.tabix_path, args.gtf_file,
                              f'chr{cromosome_num}:{tabix_query_range_start}-{tabix_query_range_end}'],
                             capture_output=True, text=True).stdout.replace("\t", "|")
    if out_str == "":
        return ""

    gene_entries_only = []
    split_out_str = out_str.split("\n")
    for line in split_out_str:  # Iterate through the lines of the output, and return the first that refers to a gene
        split_line = line.split("|")
        if len(split_line) > 2 and split_line[2] == "gene":
            gene_entries_only.append(line)

    if len(gene_entries_only) == 0:
        return ""
    elif len(gene_entries_only) == 1:
        return gene_entries_only[0]
    else:
        return find_closest_gene(lead_snp_pos, gene_entries_only)


# Only add the closest gene if the GENCODE v42 gtf file is provided. Otherwise, fill the column with NA
if args.gtf_file is not None:
    # Add a column with the whole tabix query output row
    collapsed_bed_df = collapsed_bed_df.with_columns(
        [pl.struct(["chrom", "leadSnp_pos"]).apply(get_closest_gene_tabix_output_row).alias("tabix_full_output")])

    # In those cases where there is a tabix output row, extract the gene name and add it as a column
    collapsed_bed_df = collapsed_bed_df.with_column(
        pl.when(pl.col("tabix_full_output").str.lengths() > 0)
        .then(
            pl.col("tabix_full_output").str.split("|").arr.get(-1).str.split(";").arr.get(2).str.split('"').arr.get(1))
        .otherwise("NA").alias("leadSnp_closest_gene"))

    collapsed_bed_df.drop_in_place("tabix_full_output")  # Drop the tabix output column
else:
    collapsed_bed_df = collapsed_bed_df.with_column(pl.lit("NA").alias("leadSnp_closest_gene"))

# Reorder columns
collapsed_bed_df = collapsed_bed_df.select(["chrom", "chromStart", "chromEnd", "leadSnp_pos", "leadSnp_pval",
                                            "leadSnp_closest_gene", "num_included_variants", "num_phenotypes",
                                            "phenotypes"])

########################################################################################################################
# Write the output file

output_file_name = "variant_regions.bed"
collapsed_bed_df.write_csv(args.output_filepath + '/' + output_file_name, sep="\t", has_header=True)
print(f"Done! Final regions saved to '{args.output_filepath + '/' + output_file_name}'")

# Remove the intermediate files
os.remove(args.output_filepath + "/variant_regions_step1.bed")
os.remove(args.output_filepath + "/variant_regions_step2.bed")
