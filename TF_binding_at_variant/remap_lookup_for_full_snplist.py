import os
import argparse
import polars as pl
from extract_remapdb_studies import extract_studies_for_single_snp

"""
Wrapper function to apply extract_studies_for_single_snp() to every SNP in a 'snplist' file.
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Run extract_studies_for_single_snp() for every SNP in a file")
parser.add_argument("-s", "--snplist", required=True, help="Path to file with snplist")
parser.add_argument("-r", "--remapdb", required=True, help="Path to ReMap database")
parser.add_argument("-t", "--tmpdir", required=True,
                    help="Path to temporary directory where intermediate files will be stored")
parser.add_argument("-o", "--output", required=True, help="Path to output dir. If the snplist contains several SNPs, "
                                                          "the output will be as many files, "
                                                          "named remap_studies_<chr:pos>.txt")

args = parser.parse_args()

# Validate input
if not os.path.isfile(args.snplist):
    raise ValueError("Input file does not exist")

if not os.path.isdir(args.tmpdir):
    raise ValueError("Temporary directory does not exist")

# Read snplist
df = pl.read_csv(args.snplist, sep="\t", has_header=True,
                 dtypes={"Chrom": pl.Utf8, "Pos": pl.Int64, "OA": pl.Utf8, "EA": pl.Utf8})

# Check that columns "ID", "Chrom", "Pos", "OA", "EA" exist
if df.columns != ["ID", "Chrom", "Pos", "OA", "EA"]:
    raise ValueError("Input file does not have the correct columns. The columns should be: ID, Chrom, Pos, OA, EA")

# Combine columns Chrom (without 'chr') and Pos with a ':' in between into a new column called chr_pos
df = df.with_columns([pl.format("{}:{}", pl.col("Chrom").str.replace("chr", ""), pl.col("Pos")).alias("chr_pos")])

# Run extract_studies_for_single_snp() for each chr_pos value
df = df.with_column(pl.col("chr_pos")
                    .apply(lambda s: extract_studies_for_single_snp(s,
                                                                    args.remapdb,
                                                                    args.tmpdir,
                                                                    os.path.join(args.output,
                                                                                 f"remap_studies_{str(s)}.txt"))))
