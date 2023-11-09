import argparse
import polars as pl

parser = argparse.ArgumentParser(description="Look up the FABIAN-variant output to add gain/loss information to the "
                                             "TF names in the double evidence file")

parser.add_argument("double_evidence_file", help="Path to double evidence file")
parser.add_argument("fabian_variant_output", help="Path to FABIAN-variant output file")
parser.add_argument("-m", "--map_file", required=True, help="Name of the map file between variant IDs and Fabian "
                                                            "input format.")
parser.add_argument("-o", "--output_file", required=True, help="Path to output file")

args = parser.parse_args()

# Read the double evidence file
tf_df = pl.read_csv(args.double_evidence_file, separator="\t", has_header=True)

# Split the TFs column into a list by the comma separator
tf_df = tf_df.with_columns(pl.col("TFs").str.split(",").alias("TFs"))

# Explode the TFs column
tf_df = tf_df.explode("TFs")
# Drop all nulls
tf_df = tf_df.drop_nulls(subset=["TFs"])
print(f"tf_df:\n{tf_df}")

# Drop all rows that have a missing value in the TFs column
tf_df = tf_df.drop_nulls(subset=["TFs"])

# Read the FABIAN-variant output file. This is usually a HUGE table that contains all TFs, significant or not, per SNP
fabian_df = pl.read_csv(args.fabian_variant_output, separator="\t", has_header=True)
# Fabian assigns suffixes to the variant names (i.e. .1, .2). We will remove them here
fabian_df = fabian_df.with_columns(pl.col("variant").str.split(".").map_elements(lambda s: s[0]).alias("variant"))

# Read the map file, which contains the mapping between the "variant IDs" (chr:pos) and the FABIAN input format
map_df = pl.read_csv(args.map_file, separator="\t", has_header=True, dtypes={"ID": pl.Utf8, "Chrom:PosOA>EA": pl.Utf8})
print(f"map_df:\n{map_df}")

# Add an ID column to fabian_df by matching the Chrom:PosOA>EA column
fabian_df = fabian_df.join(map_df, left_on="variant", right_on="Chrom:PosOA>EA", how="inner")
# Rename the 'tf' column to 'TFs' to match the double evidence file
fabian_df = fabian_df.rename({"tf": "TFs"})
print(f"fabian_df:\n{fabian_df}")

# Extract final table. Join by ID and TF
out_df = tf_df.join(fabian_df, on=["ID", "TFs"], how="inner")
# Keep only the columns ID, variant, TFs, prediction and score
out_df = out_df.select(["ID", "variant", "TFs", "prediction", "score"])
# Add column 'abs_score' that contains the absolute value of the 'score' column
out_df = out_df.with_columns(pl.col("score").abs().alias("abs_score"))
# For the rows that are identical except for the 'score' column, keep the one with the highest absolute score
out_df = out_df.sort(by="abs_score", descending=True).unique(subset=["ID", "variant", "TFs"], keep="first")
# Drop abs_score column
out_df = out_df.drop("abs_score")

# Use regex to extract the 'letters>letters' part from the 'variant' column
out_df = out_df.with_columns(pl.col("variant").str.extract(r"([ATCG]+>[ATCG])").alias("OA>EA"))
# Separate OA and EA into two columns
out_df = out_df.with_columns([pl.col("OA>EA").str.split(by=">").map_elements(lambda s: s[0]).alias("OA"),
                              pl.col("OA>EA").str.split(by=">").map_elements(lambda s: s[1]).alias("EA")])
# Drop the 'OA>EA' and the 'variant' columns
out_df = out_df.drop(["OA>EA", "variant"])

# FINAL SORTING - we will extract chr and pos from ID, and sort by chr, pos and score
# Extract chr and pos from ID
out_df = out_df.with_columns([pl.col("ID").str.split(by=":").map_elements(lambda s: s[0]).str.lstrip("chr").cast(pl.Int32).alias("Chrom"),
                              pl.col("ID").str.split(by=":").map_elements(lambda s: s[1]).cast(pl.Int32).alias("Pos")])
# Sort
out_df = out_df.sort(by=["Chrom", "Pos", "score"], descending=[False, False, True])

# Change the final order of the columns to: ID, OA, EA, TFs, prediction, score
out_df = out_df.select(["ID", "OA", "EA", "TFs", "prediction", "score"])
out_df.write_csv(args.output_file, separator="\t", has_header=True)

print(f"out_df:\n{out_df}")