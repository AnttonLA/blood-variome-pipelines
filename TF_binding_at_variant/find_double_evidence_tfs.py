import argparse
import sys
import os
import polars as pl

"""
This script will look at the output of both the ReMap database lookup and the predicted transcription factor binding
motif disruption (FABIAN-Variant), and it will find the TFs that show up in both.
There is a good chance that this TFs are involved in the effect of the variant, since there is both ChIP-Seq
and motif disruption evidence for them.
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Look at ReMap and TFBS disruption prediction (FABIAN-Variant) and find "
                                             "TFs that show up in both")
parser.add_argument("-r", "--remap_file", required=True, help="Path to directory ReMap output files are stored")
parser.add_argument("-f", "--fabian_table_file", required=True, help="Path to fabian output file")
parser.add_argument("-o", "--output_file", required=True, help="Output file")

args = parser.parse_args()

# Input sanitation. Make sure input files exist.
if not os.path.isfile(args.remap_file) or os.path.getsize(args.remap_file) == 0:
    raise ValueError(f"ReMap file {args.remap_file} does not exist or is empty")

elif args.fabian_table_file and not os.path.isfile(args.fabian_table_file):
    raise ValueError(f"Fabian table file {args.fabian_table_file} does not exist")


def print_status(percent):
    """Prints a status bar to the console. Used to show progress of the scrip."""
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()


# Read in ReMap file
remap_df = pl.read_csv(args.remap_file, separator="\t", has_header=True)
remap_df = remap_df.select(["chr", "pos", "transcription_factor"])
remap_df = remap_df.rename({"chr": "Chr", "pos": "Pos", "transcription_factor": "transcription_factor"})
print(f"remap_df:\n{remap_df}")

# Read in FABIAN-Variant output file
fabian_df = pl.read_csv(args.fabian_table_file, separator="\t", has_header=True)
fabian_df = fabian_df.rename({"Chrom": "Chr", "Pos": "Pos", "TF": "transcription_factor"})
print(f"fabian_df:\n{fabian_df}")

# TODO: all of this renaming could be avoided if we do it during the table creation...

df = fabian_df.join(remap_df, on=["Chr", "Pos", "transcription_factor"], how="inner")
# Remove duplicate rows
df = df.unique()

# If 'effect' is not already a column, add it. 'gain' if the score is positive and 'loss' if the score is negative
if "effect" not in df.columns:
    df = df.with_columns(pl.when(pl.col("score") > 0).then(pl.lit("gain")).otherwise(pl.lit("loss")).alias("effect"))
# Sort by Chr and Pos
df = df.sort(by=["Chr", "Pos"])

print(f"df:\n{df}")
df.write_csv(args.output_file, separator="\t", has_header=True)
