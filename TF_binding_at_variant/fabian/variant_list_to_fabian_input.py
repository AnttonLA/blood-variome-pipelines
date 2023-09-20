import polars as pl

# Read the input TSV file into a Polars DataFrame
input_file = "/home/antton/Projects/blood-variome-pipelines/TF_binding_at_variant/tests/test_data/dummy_variant_list.tsv"
df = pl.read_csv(input_file, sep='\t')

# Create a new column with the desired format
df = df.with_column((pl.col("Chrom") + ":" + pl.col("Pos") + pl.col("OA") + ">" + pl.col("EA")).alias("Chrom_Pos_OA_EA"))

# Select only the new column and convert it to a list
output_series = df.select(["Chrom_Pos_OA_EA"])

# Write the list to a file
output_file = "/home/antton/Projects/blood-variome-pipelines/TF_binding_at_variant/tests/test_data/dummy_fabian_input.txt"

output_series.write_csv(output_file, sep='\t', has_header=False)
