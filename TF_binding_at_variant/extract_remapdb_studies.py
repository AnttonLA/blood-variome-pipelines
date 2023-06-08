import polars as pl
import os
import argparse
import sys

"""
This script will take the position of a SNP ('snp_full_pos') and produce a tabix query to the filepath stored in
'remap_file'. It will then process the output of the tabix query to produce a table with the following columns:
- study_accession
- transcription_factor
- biotype
- distance_to_snp

"""

# Parse arguments
parser = argparse.ArgumentParser(description='Take the position of a SNP and produce a tabix query to the ReMap '
                                             'dabase file.')
parser.add_argument('snp_pos', type=str, help='Position of SNP, e.g. 12:62910804')
parser.add_argument('-r', '--remap_file', type=str, help='Path to the ReMap database file', required=True)
parser.add_argument('-o', '--output', type=str, help='Name of the output file', required=True)

args = parser.parse_args()

# Input sanitation
if not ':' in args.snp_pos:
    raise ValueError("SNP position must be in the format '<chr>:<position>'")

if not args.snp_pos.split(":")[1].isdigit():
    raise ValueError("SNP position must be in the format '<chr>:<position>'")

if not args.output.endswith(".csv"):
    sys.stderr.write("WARNING: Output file must be a CSV file. Adding extension '.csv' to output file name.\n")
    args.output = args.output + ".csv"

# Create 'snp_full_pos' variable
chr = args.snp_pos.split(":")[0]  # Chromosome of the SNP
pos = args.snp_pos.split(":")[1]  # Position of the SNP
snp_full_pos = f"chr{chr}:{pos}-{pos}"  # Full position of the SNP in format of chr12:62910804-62910804, for tabix query

print(f"Requested position: <{snp_full_pos}>")

# Check if ReMap file exists
if not os.path.isfile(args.remap_file):
    raise ValueError("ReMap file does not exist")

# Check it is bgzipped
if not args.remap_file.endswith(".gz"):
    raise ValueError("ReMap file must be bgzipped")

# Check that the tabix index file exists
if not os.path.isfile(args.remap_file + ".tbi"):
    raise ValueError("No tabix index file could be found for the ReMap file.")

remap_file = args.remap_file  # Full ReMap BED file, downloaded from the ReMap website

# Compose tabix query
command = f"tabix {remap_file} {snp_full_pos}"

# Run command and store output. Requires tabix.
output = os.popen(command).read()

# Check if the output is empty
if output == "":
    raise ValueError("No ReMap entries found for this SNP.")

# Save the output to a file
tmp_file = f"data/tabix_slices/remap2022_all_macs2_hg38_v1_0.{snp_full_pos}.bed"
with open(tmp_file, "w") as f:
    f.write(output)

df = pl.read_csv(tmp_file, sep="\t",
                 new_columns=["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"])

# Distance of SNP to peak. We'll use this to sort entries later.
# It'd probably be fine to calculate distance to thickStart, but just to be safe we'll calculate the center of the peak.
df = df.with_column(pl.struct(['thickStart', 'thickEnd'])
                    .apply(lambda s: int((s['thickStart'] + s['thickEnd']) / 2))
                    .alias("thickCenter"))

# Calculate distance to SNP pos (stored in variable 'pos')
df = df.with_column(pl.col("thickCenter").apply(lambda s: abs(int(s) - int(pos))).alias("distance_to_peak"))

# Split the name column by '.' and make new columns out of the three entries
df = df.with_columns([pl.col("name").str.split(".").apply(lambda s: s[0]).alias("study_accession"),
                      pl.col("name").str.split(".").apply(lambda s: s[1]).alias("transcription_factor"),
                      pl.col("name").str.split(".").apply(lambda s: s[2]).alias("biotype")])

out_df = df.select(["study_accession", "transcription_factor", "biotype", "distance_to_peak"])

# Sort by distance from SNP to ChIP-seq peak.
out_df = out_df.sort("distance_to_peak")

# Write to file
out_df.write_csv(args.output)
