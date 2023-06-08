import argparse
import polars as pl
import os

"""
This script will look at the output of both the ReMap database lookup and the perfectos-ape run and find the TFs that
show up in both.
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Look at ReMap and perfectos-ape output and find TFs that show up in both")
parser.add_argument("--remap", required=True, help="(REQUIRED) Path to ReMap output file")
parser.add_argument("--perfectos", required=True, help="(REQUIRED) Path to perfectos-ape output file")

args = parser.parse_args()

# Input sanitation. Make sure input files exist.
if not os.path.isfile(args.remap):
    raise ValueError("ReMap file does not exist")

if not os.path.isfile(args.perfectos):
    raise ValueError("perfectos-ape file does not exist")

# Read in ReMap file
remap_df = pl.read_csv(args.remap, sep=",")

# Read in perfectos-ape file
perfectos_df = pl.read_csv(args.perfectos, sep="\t", has_header=True)

# Get all the unique values of the second column of remap_df
remap_tfs = remap_df.select("transcription_factor").unique().get_column("transcription_factor").to_list()

print(f"ReMap TFs: {remap_tfs}")

# Split all the values of the second column of perfectos_df by _ and take the first element of each split
perfectos_df = perfectos_df.with_column(pl.col("motif").str.split("_").apply(lambda s: s[0]).alias("tf_only"))
perfectos_tfs = perfectos_df.select("tf_only").unique().get_column("tf_only").to_list()

print(f"Perfectos TFs: {perfectos_tfs}")

# Find the intersection of the two lists
output_tfs = list(set(remap_tfs).intersection(perfectos_tfs))

print("The TFs that show up in both ReMap and Hocomoco are: ", output_tfs)
