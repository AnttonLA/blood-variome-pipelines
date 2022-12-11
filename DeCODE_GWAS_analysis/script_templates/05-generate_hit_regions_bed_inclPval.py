import polars as pl
import os
import subprocess

"""
This script is used to generate a bed file with the regions of the hits table that are in within 1Mb of one other.
It also includes the phenotypes of the hits in those regions.

Note that it uses 'bedtools'. You will need to edit the path in the subprocess call to wherever your bedtools is.
"""

folder_path = "/media/nvme/antton/Frequency_and_Ratio/combined_output/all_significant_variant_tables/"

# Create empty output df, specify datatypes of each column.
bed_df = pl.DataFrame({"chrom": [], "chromStart": [], "chromEnd": [], "lowest_pval": [], "phenotypes": []},
                      columns=[("chrom", pl.Int64), ("chromStart", pl.Int64), ("chromEnd", pl.Int64),
                               ("lowest_pval", pl.Float64), ("phenotypes", str)])

print("Welcome!\nSTEP 1: Reading files in the 'all_significant_variant_tables' folder.")
for i, file in enumerate(os.listdir(folder_path)):
    df = pl.read_csv(folder_path + file, sep="\t", columns=["pval", "chromosome", "position", "phenotype"])
    # Replace "chrX" with "chr23"
    df = df.with_column(
        pl.when(pl.col("chromosome") == "chrX").then("chr23").otherwise(pl.col("chromosome")).alias("chromosome"))
    for row in df.rows():
        # ("pval", "chromosome", "position", "phenotype")
        if bed_df.is_empty():  # First entry
            bed_df = pl.DataFrame({"chrom": [int(row[1][3:])],
                                   "chromStart": [max(2, row[2] - 10 ** 6)],
                                   "chromEnd": [row[2] + 10 ** 6],
                                   "pval": [row[0]],
                                   "phenotypes": [row[3]]})
        else:
            bed_df = bed_df.extend(pl.DataFrame({"chrom": [int(row[1][3:])],
                                                 "chromStart": [max(2, row[2] - 10 ** 6)],
                                                 "chromEnd": [row[2] + 10 ** 6],
                                                 "pval": [row[0]],
                                                 "phenotypes": [row[3]]}))
    if not (i % 25):
        percent = 100 * (i / len(os.listdir(folder_path)))
        print("#" * int(percent), f"({round(percent, 2)}%)")
    # if i >= 7:
    #    break

bed_df = bed_df.sort([pl.col("chrom"), pl.col("chromStart")])
bed_df.write_csv("variant_regions_step1.bed", sep="\t", has_header=True)

print("Completed STEP 1!\nSTEP 2: Merging all overlapping regions with 'bedtools'.")
# Use bedtools to collapse the overlapping ranges while keeping the phenotype names
with open("variant_regions_step2.bed", "w") as f:
    subprocess.run(["/home/antton/Programs/bedtools2/bin/bedtools", "merge", "-i", "variant_regions_step1.bed",
                    "-header", "-c", "4,5", "-o", "distinct"],
                   check=True, stdout=f)

# Read the collapsed bed file and convert it to a polars dataframe
collapsed_bed_df = pl.read_csv("variant_regions_step2.bed", sep="\t", columns=["chrom", "chromStart", "chromEnd",
                                                                               "pval", "phenotypes"])

# Add a column with the numbers of phenotypes with a significant p-value within that region.
collapsed_bed_df.insert_at_idx(3, collapsed_bed_df.select("phenotypes").apply(
    lambda x: len(str(x).split(',')) - 1).to_series())
collapsed_bed_df = collapsed_bed_df.rename({"apply": "num_phenotypes"})

print(collapsed_bed_df["pval"][0])
# Add a column with the lowest pvalue for a variant in the region
collapsed_bed_df.replace("pval", collapsed_bed_df.select("pval").apply(
    lambda x: "{:.2e}".format(min([float(i) for i in str(x[-1]).split(',')]))).to_series())
collapsed_bed_df = collapsed_bed_df.rename({"pval": "lowest_pval"})

output_file_name = "variant_regions_wPvals.bed"
collapsed_bed_df.write_csv(output_file_name, sep="\t", has_header=True)
print(f"Done! Final regions saved to '{output_file_name}'")
