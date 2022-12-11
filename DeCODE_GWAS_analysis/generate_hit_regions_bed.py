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
bed_df = pl.DataFrame({"chrom": [], "chromStart": [], "chromEnd": [], "lowest_pval": [], "phenotypes": []},
                      columns=[("chrom", pl.Int64), ("chromStart", pl.Int64), ("chromEnd", pl.Int64),
                               ("lowest_pval", pl.Float64), ("phenotypes", str)])

print("Welcome!\nReading all summary stats files in the provided folder...")
for i, file in enumerate(os.listdir(args.folder_path)):
    df = pl.read_csv(args.folder_path + file, sep="\t", columns=["pval", "chromosome", "position", "phenotype"])
    for row in df.rows():
        # ("pval", "chromosome", "position", "phenotype")
        if bed_df.is_empty():  # First entry
            bed_df = pl.DataFrame({"chrom": [int(row[1])],
                                   "chromStart": [max(2, row[2] - 10 ** 6)],
                                   "chromEnd": [row[2] + 10 ** 6],
                                   "pval": [row[0]],
                                   "phenotypes": [row[3]]})
        else:
            bed_df = bed_df.extend(pl.DataFrame({"chrom": [int(row[1])],
                                                 "chromStart": [max(2, row[2] - 10 ** 6)],
                                                 "chromEnd": [row[2] + 10 ** 6],
                                                 "pval": [row[0]],
                                                 "phenotypes": [row[3]]}))
    if not (i % 10):
        percentage = 100 * (i / len(os.listdir(args.folder_path)))
        print_status(percentage)
    # if i >= 7:
    #    break

bed_df = bed_df.sort([pl.col("chrom"), pl.col("chromStart")])
bed_df.write_csv(args.output_filepath + "/variant_regions_step1.bed", sep="\t", has_header=True)

print("Complete!\nMerging all overlapping regions with 'bedtools'.")
# Use bedtools to collapse the overlapping ranges while keeping the phenotype names
with open(args.output_filepath + "/variant_regions_step2.bed", "w") as f:
    subprocess.run(["/home/antton/Programs/bedtools2/bin/bedtools", "merge", "-i",
                    f"{args.output_filepath}/variant_regions_step1.bed", "-header", "-c", "4,5", "-o", "distinct"],
                   check=True, stdout=f)

# Read the collapsed bed file and convert it to a polars dataframe
collapsed_bed_df = pl.read_csv(args.output_filepath + "/variant_regions_step2.bed", sep="\t",
                               columns=["chrom", "chromStart", "chromEnd", "pval", "phenotypes"])

# Add a column with the numbers of phenotypes with a significant p-value within that region.
collapsed_bed_df.insert_at_idx(3, collapsed_bed_df.select("phenotypes").apply(
    lambda x: len(str(x).split(',')) - 1).to_series())
collapsed_bed_df = collapsed_bed_df.rename({"apply": "num_phenotypes"})

# Add a column with the lowest pvalue for a variant in the region
collapsed_bed_df.replace("pval", collapsed_bed_df.select("pval").apply(
    lambda x: "{:.2e}".format(min([float(i) for i in str(x[-1]).split(',')]))).to_series())
collapsed_bed_df = collapsed_bed_df.rename({"pval": "lowest_pval"})

output_file_name = "variant_regions.bed"
collapsed_bed_df.write_csv(args.output_filepath + '/' + output_file_name, sep="\t", has_header=True)
print(f"Done! Final regions saved to '{args.output_filepath + '/' + output_file_name}'")

# Remove the intermediate files
os.remove(args.output_filepath + "/variant_regions_step1.bed")
os.remove(args.output_filepath + "/variant_regions_step2.bed")
