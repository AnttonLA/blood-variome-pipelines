import argparse
import sys

import polars as pl
import os

"""
This script will look at the output of both the ReMap database lookup and the predicted transcription factor binding
motif disruption (either FABIAN-Variant or PERFECTOS-APE), and it will find the TFs that show up in both.
There is a good chance that this TFs are involved in the effect of the variant, since there is both Chip-Seq
and motif disruption evidence for them.
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Look at ReMap and TFBS disruption prediction (either PERFECTOS-APE or "
                                             "FABIAN-Variant) and find TFs that show up in both")
parser.add_argument("-r", "--remap_dir", required=True, help="Path to directory ReMap output files are stored")
tf_disruption_group = parser.add_mutually_exclusive_group(required=True)
tf_disruption_group.add_argument("-p", "--perfectos", help="Path to perfectos-ape output file")
tf_disruption_group.add_argument("-f", "--fabian", help="Path to fabian output file")
parser.add_argument("-fm", "--fabian_map", required=False, help="Path to map file between variant IDs and Fabian "
                                                                "input in chr:posOA>EA (cpoa) format. Only needed if "
                                                                "Fabian was chosen for the TFBS disruption prediction.")
parser.add_argument("-ft", "--fabian_threshold", required=False, default=0.2, help="Threshold for Fabian score. Value "
                                                                                   "between 0 and 1. Default: 0.2")
parser.add_argument("-o", "--output", required=True, help="Output file")

args = parser.parse_args()

# Input sanitation. Make sure input files exist.
if not os.path.isdir(args.remap_dir):
    raise ValueError(f"ReMap dir {args.remap_dir} does not exist")
# If it does not end in a '/', add it
if args.remap_dir[-1] != '/':
    args.remap_dir += '/'

if args.perfectos and not os.path.isfile(args.perfectos):
    raise ValueError(f"perfectos-ape file {args.perfectos} does not exist")
elif args.fabian and not os.path.isfile(args.fabian):
    raise ValueError(f"Fabian file {args.fabian} does not exist")

# If fabian was chosen, make sure a map file was provided
if args.fabian and not args.fabian_map:
    raise ValueError("Fabian was chosen for the TFBS disruption prediction, but no map file was provided")

# If fabian was chosen, check the map file exists
if args.fabian and not os.path.isfile(args.fabian_map):
    raise ValueError(f"Fabian map file {args.fabian_map} does not exist")

# If perfectos-ape was but fabian map or threshold were provided, raise an warning that they will be ignored
if args.perfectos and (args.fabian_map or args.fabian_threshold != 0.2):
    print("Warning: perfectos-ape was chosen as the TFBS disruption prediction method, but a Fabian map file or "
          "threshold were provided. These will be ignored.")


def print_status(percent):
    """Prints a status bar to the console. Used to show progress of the script when we iterate through files."""
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()


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
sys.stdout.write("\nLoading ReMap data...\n")
remap_dict = {}
variant_list = remap_df.select("ID").unique().get_column("ID").to_list()
print("REMAP DATAFRAME\n", remap_df)
print("VARIANT LIST\n:", variant_list[:2])
for i, variant in enumerate(variant_list):
    # Progress bar
    percentage = int(i / len(variant_list) * 100)
    print_status(percentage)

    remap_dict[variant] = remap_df.filter(pl.col("ID") == variant).select("transcription_factor") \
        .get_column("transcription_factor").to_list()


# ReMap dict IDs are chr:pos only

def process_perfectos_output() -> dict:
    """
    Read in the output of perfectos-ape and return a dict where the keys are variants and the values are lists of TFs
    that correspond to that variant.

    :return: tf_disruption_dict
    """
    # Read in perfectos-ape file
    perfectos_df = pl.read_csv(args.perfectos, sep="\t", has_header=True)

    # Split all the values of the second column of perfectos_df by _ and take the first element of each split (TF name)
    perfectos_df = perfectos_df.with_column(pl.col("motif").str.split("_").apply(lambda s: s[0]).alias("tf_only"))

    # Add 'ID' column to perfectos_df for mapping between <chr:pos> and <chr:pos_OA_EA>
    perfectos_df = perfectos_df.with_column(pl.col("# SNP name").str.split("_").apply(lambda s: s[0]).alias("ID"))

    # Create dict where variants are keys and a list of TFs is the value, same as we did for the ReMap dict
    perfectos_dict = {}
    variant_list = perfectos_df.select("# SNP name").unique().get_column("# SNP name").to_list()

    for variant in variant_list:
        perfectos_dict[variant] = perfectos_df.filter(pl.col("# SNP name") == variant).select("tf_only") \
            .get_column("tf_only").to_list()

    return perfectos_dict


def process_fabian_output() -> dict:
    """
    Read in the output of Fabian and return a dict where the keys are variants and the values are lists of TFs
    that correspond to that variant.

    :return: fabian_dict - {'chr:pos_OA_EA': [TF1, TF2, ...], ...}
    """
    sys.stdout.write("\nCommencing Fabian data processing...\n")

    # Read in Fabian file
    fabian_df = pl.read_csv(args.fabian, sep="\t", has_header=True)

    # Filter out rows where the value of teh "prediction" column is "NA"
    fabian_df = fabian_df.filter(pl.col("prediction") != "NA")

    # Filter out rows where the absolute value of the "score" column is less than a threshold
    fabian_score_threshold = args.fabian_threshold  # TODO: this threshold is arbitrary. Some justification is needed.
    fabian_df = fabian_df.filter(pl.col("score").abs() >= fabian_score_threshold)

    # Fabian assigns suffixes to the variant names (i.e. .1, .2). We will remove them here
    fabian_df = fabian_df.with_column(pl.col("variant").str.split(".").apply(lambda s: s[0]).alias("variant"))
    # We will now replace the keys of the dict with the same format as the ReMap dict. Use the map file to do this.
    fabian_format_to_id_map_df = pl.read_csv(args.fabian_map, sep='\t', has_header=True).select(
        ["ID", "Chrom:PosOA>EA"])
    fabian_format_to_id_map = dict(
        zip(fabian_format_to_id_map_df.select("Chrom:PosOA>EA").get_column("Chrom:PosOA>EA").to_list(),
            fabian_format_to_id_map_df.select("ID").get_column("ID").to_list()))
    # Drop any rows where the variant is not in the map file
    fabian_df = fabian_df.filter(pl.col("variant").is_in(pl.Series('variant', list(fabian_format_to_id_map.keys()))))
    # Replace the keys
    fabian_df = fabian_df.with_column(pl.col("variant").apply(lambda s: fabian_format_to_id_map[s]).alias("variant"))

    fabian_dict = {}
    # Separate the df into sub-dfs based on the unique values of the "variant" column
    variant_list = fabian_df.select("variant").unique().get_column("variant").to_list()
    for i, variant in enumerate(variant_list):
        # Progress bar
        percent = int(i / len(variant_list) * 100)
        print_status(percent)
        variant_short = variant.split("_")[
            0]  # Remove the _<OA>_<EA> suffix from the dict key (i.e. chr1:23935190_G_T to chr1:23935190)
        # Make a list out of every unique entry of the 'tf' column
        fabian_dict[variant_short] = fabian_df.filter(pl.col("variant") == variant).select("tf").unique().get_column(
            "tf").to_list()

    return fabian_dict


if args.perfectos:
    tf_disruption_dict = process_perfectos_output()
elif args.fabian:
    tf_disruption_dict = process_fabian_output()
else:
    raise ValueError("Neither perfectos-ape nor fabian were chosen as the TFBS disruption prediction method. This "
                     "should never happen.")

print(f"TF disruption dict: \n{list(tf_disruption_dict.items())[:2]}")
print(f"Remap dict: \n{list(remap_dict.items())[:2]}")
# Find the intersection of the two dicts
sys.stdout.write("\nPerforming cross-check...\n")

output_dict = {}
for i, variant in enumerate(remap_dict.keys()):
    percentage = int(i / len(remap_dict.keys()) * 100)
    print_status(percentage)
    if variant in list(tf_disruption_dict.keys()):
        if 'X' in variant:
            variant_sanitized = variant.replace("X", "23")
        else:
            variant_sanitized = variant
        output_dict[variant_sanitized] = list(set(remap_dict[variant]).intersection(set(tf_disruption_dict[variant])))


# Define a custom sorting key function to deal with genomic coordinate sorting
def genome_sort_key(position):
    """
    Tiny util to sort genomic coordinates properly
    :param position:
    :return:
    """
    # Split the position into chromosome and position parts
    parts = position.split(':')
    chromosome = int(parts[0].lstrip('chr'))  # Remove 'chr' and convert to int
    position = int(parts[1])
    return chromosome, position


# Print the first 10 elements of 'output_dict'
print("\nOutput_dict:\n", list(output_dict.items())[:4])

# Write output
with open(args.output, "w") as f:
    f.write("ID\tTFs\n")
    if list(output_dict.keys())[0].startswith("chr"):
        for variant in sorted(output_dict.keys(), key=lambda x: genome_sort_key(x)):
            f.write(variant + "\t" + ",".join(output_dict[variant]) + "\n")
    else:  # Sort normally if the keys are not in chr:pos format (i.e. they are proper IDs)
        for variant in sorted(output_dict.keys()):
            f.write(variant + "\t" + ",".join(output_dict[variant]) + "\n")
