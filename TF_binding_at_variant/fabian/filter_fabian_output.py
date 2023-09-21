import polars as pl

in_file = "/home/antton/Projects/blood-variome-pipelines/TF_binding_at_variant/tests/test_data/dummy_fabian_output.tsv"

df = pl.read_csv(in_file, sep='\t', has_header=True)

# Filter out rows where the value of teh "prediction" column is "NA"
df = df.filter(pl.col("prediction") != "NA")

# Filter out rows where the absolute value of the "score" column is less than 0.1 # TODO: this threshold is arbitrary
df = df.filter(pl.col("score").abs() >= 0.1)

fabian_dict = {}
# Separate the df into sub-dfs based on the unique values of the "variant" column
variant_list = df.select("variant").unique().get_column("variant").to_list()
for variant in variant_list:
    # Make a list out of every unique entry of the 'tf' column
    fabian_dict[variant] = df.filter(pl.col("variant") == variant).select("tf").unique().get_column("tf").to_list()

print(fabian_dict)
