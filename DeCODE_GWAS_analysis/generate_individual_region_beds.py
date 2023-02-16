import os
import polars as pl
import argparse

"""
This script will go through a BED-like file specifying the regions in the genome where we have GWAS hits, and will check
all of the GWAS hits in the hits_only_sumstats folder for hits within that region. It will then output a BED file for
each of the regions with positions, effect alleles, p-values and phenotypes of the hits within that region.
"""

parser = argparse.ArgumentParser(description="Make a BED file for each of the regions specified in the"
                                             "'variant_regions.bed' file. Each file will contain the positions, effect"
                                             "alleles, p-values and phenotypes of the hits within that region.")
parser.add_argument("--regions_file", type=str, required=True,
                    help="Path to the BED-like file containing the regions in the genome where we have GWAS hits.")
parser.add_argument("--hits_bed_files_folder", type=str, required=True,
                    help="Path to the folder containing the GWAS hit tables.")
parser.add_argument("--output_folder", type=str, required=True,
                    help="Path to the folder where the output files will be written to.")

args = parser.parse_args()

# Check if the paths are valid
if not os.path.exists(args.regions_file):
    raise ValueError("The given regions file does not exist.")
if not os.path.exists(args.hits_bed_files_folder):
    raise ValueError("The given hits BED files folder does not exist.")
# Check if the output folder exists, if not create it
if not os.path.exists(args.output_folder):
    os.mkdir(args.output_folder)

# If folder names do not end with a slash, add it
if args.hits_bed_files_folder[-1] != "/":
    args.hits_bed_files_folder += "/"
if args.output_folder[-1] != "/":
    args.output_folder += "/"

# Load regions file
df = pl.read_csv(args.regions_file, sep='\t', columns=["chrom", "chromStart", "chromEnd", "phenotypes"])
df = df.with_column(pl.col("phenotypes").str.split(by=","))  # Split phenotypes column into a list

# Iterate row by row through the regions file
for region_number in range(len(df)):
    region_chrom = df[region_number, 0]
    region_start = df[region_number, 1]
    region_end = df[region_number, 2]
    pheno_list = df[region_number, 3].to_list()  # List of phenotypes in the selected region

    final_df = pl.DataFrame(columns=[("chromosome", pl.Int32),
                                     ("position", pl.Int32),
                                     ("EA", pl.Utf8),
                                     ("pval", pl.Float64),
                                     ("phenotype", pl.Utf8)])
    # Read, filter and stack the "hit_tables" of each phenotype taking only hits in the selected region
    for pheno in pheno_list:
        hits_bed_file = args.hits_bed_files_folder + pheno + ".txt"
        hits_bed_df = pl.read_csv(hits_bed_file, sep='\t',
                                  columns=["pval", "chromosome", "position", "EA", "phenotype"],
                                  dtypes={"chromosome": pl.Int32, "position": pl.Int32, "EA": pl.Utf8,
                                          "pval": pl.Float64, "phenotype": pl.Utf8})
        # Reorder columns
        hits_bed_df = hits_bed_df.select(["chromosome", "position", "EA", "pval", "phenotype"])
        hits_bed_df = hits_bed_df.filter(hits_bed_df["chromosome"] == df[region_number, 0])  # Same chromosome
        hits_bed_df = hits_bed_df.filter(hits_bed_df["position"] >= df[region_number, 1])  # Within region
        hits_bed_df = hits_bed_df.filter(hits_bed_df["position"] <= df[region_number, 2])

        final_df = pl.concat([final_df, hits_bed_df])
    final_df.write_csv(args.output_folder + f"region_{region_number+1}_chr{region_chrom}:{region_start}-{region_end}.bed",
                       sep='\t', has_header=True)
