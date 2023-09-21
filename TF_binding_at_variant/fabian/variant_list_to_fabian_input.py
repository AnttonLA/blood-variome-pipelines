import os
import polars as pl
import argparse

parser = argparse.ArgumentParser(description="Convert a variant list to the format required by FABIAN-Variant")
parser.add_argument("input_file", help="Path to input variant file")
parser.add_argument("-m", "--map_file", required=False, help="Name of the map file between variant IDs and Fabian "
                                                             "input format. If none is provided, the output file name"
                                                             "will be used with the extension .map")
parser.add_argument("-o", "--output_file", required=True, help="Path to output file")

args = parser.parse_args()

# Check if the input file exists
if not os.path.isfile(args.input_file):
    raise ValueError(f"Input file {args.input_file} does not exist")

# Ensure the input file has the necessary columns (ID, Chrom, Pos, OA, EA)
df = pl.read_csv(args.input_file, sep='\t', has_header=True)
if not {"ID", "Chrom", "Pos", "OA", "EA"}.issubset(set(df.columns)):
    raise ValueError(f"Input file {args.input_file} does not have the necessary columns (ID, Chrom, Pos, OA, EA)")

# Create a new column with the desired format
df = df.with_column((pl.col("Chrom") + ":" + pl.col("Pos") + pl.col("OA") + ">" + pl.col("EA")).alias("Chrom_Pos_OA_EA"))

# Produce a map where the keys are "ID" and the values are "Chrom_Pos_OA_EA". Save to file.
df = df.select(["ID", "Chrom_Pos_OA_EA"])
map_file_name = args.map_file if args.map_file else args.output_file + ".map"
df.write_csv(map_file_name, sep='\t', has_header=True)

# Select only the new column and convert it to a list
output_series = df.select("Chrom_Pos_OA_EA")

# Write the list to a file
output_series.write_csv(args.output_file, sep='\t', has_header=False)
