import argparse
import polars as pl
import os

"""
This script will look at the output of both the ReMap database lookup and the perfectos-ape run and find the TFs that
show up in both.
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Look at ReMap and perfectos-ape output and find TFs that show up in both")
parser.add_argument("-r", "--remap_dir", required=True, help="Path to directory ReMap output files are stored")
parser.add_argument("-p", "--perfectos", required=True, help="Path to perfectos-ape output file")
parser.add_argument("-o", "--output", required=True, help="Output file")

args = parser.parse_args()

# Input sanitation. Make sure input files exist.
if not os.path.isdir(args.remap_dir):
    raise ValueError(f"ReMap dir {args.remap_dir} does not exist")

if not os.path.isfile(args.perfectos):
    raise ValueError(f"perfectos-ape file {args.perfectos} does not exist")

# Read in ReMap files
remap_df = pl.DataFrame(columns={"study_accession": pl.Utf8,
                                 "transcription_factor": pl.Utf8,
                                 "biotype": pl.Utf8,
                                 "distance_to_peak": pl.Int64,
                                 "ID": pl.Utf8})

for file in os.listdir(args.remap_dir):
    if file.startswith("remap_studies_") and file.endswith(".txt"):
        variant_chr_pos = file.replace("remap_studies_", "").replace(".txt", "")
        variant_id = "chr" + variant_chr_pos.split(":")[0] + ":" + variant_chr_pos.split(":")[1]  # chr<chr>:<pos>
        sub_remap_df = pl.read_csv(args.remap_dir + file, sep="\t", dtypes={"distance_to_peak": pl.Int64})
        sub_remap_df = sub_remap_df.with_column(pl.lit(variant_id).alias("ID"))
        remap_df = remap_df.vstack(sub_remap_df)

remap_df = remap_df.select(["transcription_factor", "ID"])

# Create a dict where the keys are unique IDs and the values are lists of TFs that correspond to that ID
remap_dict = {}
variant_list = remap_df.select("ID").unique().get_column("ID").to_list()
for variant in variant_list:
    remap_dict[variant] = remap_df.filter(pl.col("ID") == variant).select("transcription_factor") \
        .get_column("transcription_factor").to_list()

# Read in perfectos-ape file
perfectos_df = pl.read_csv(args.perfectos, sep="\t", has_header=True)

# Split all the values of the second column of perfectos_df by _ and take the first element of each split (TF name)
perfectos_df = perfectos_df.with_column(pl.col("motif").str.split("_").apply(lambda s: s[0]).alias("tf_only"))

# Add 'ID' column to perfectos_df for mapping between <chr:pos> and <chr:pos_OA_EA>
perfectos_df = perfectos_df.with_column(pl.col("# SNP name").str.split("_").apply(lambda s: s[0]).alias("ID"))

# Create the map out of two lists
chr_pos_list = perfectos_df.select("ID").get_column("ID").to_list()
chr_pos_oa_ea_list = perfectos_df.select("# SNP name").get_column("# SNP name").to_list()
chr_pos_oa_ea_map = {chr_pos_list[i]: chr_pos_oa_ea_list[i] for i in range(len(chr_pos_list))}

# Create dict where variants are keys and a list of TFs is the value, same as we did for the ReMap dict
perfectos_dict = {}
variant_list = perfectos_df.select("# SNP name").unique().get_column("# SNP name").to_list()

for variant in variant_list:
    variant_no_alleles = variant.split("_")[0]  # Remove _EA_OA from the end of the variant name so that names match the ReMao dict
    perfectos_dict[variant_no_alleles] = perfectos_df.filter(pl.col("# SNP name") == variant).select("tf_only") \
        .get_column("tf_only").to_list()


# Find the intersection of the two dicts
output_dict = {}
for variant in remap_dict.keys():
    if variant in perfectos_dict.keys():
        output_dict[variant] = list(set(remap_dict[variant]).intersection(set(perfectos_dict[variant])))

# Write output
with open(args.output, "w") as f:
    f.write("ID\tTFs\n")
    for variant in output_dict.keys():
        f.write(chr_pos_oa_ea_map[variant] + "\t" + ",".join(output_dict[variant]) + "\n")
