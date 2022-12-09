import os

import polars as pl
import argparse

parser = argparse.ArgumentParser(description="Take a template summary statistics file and swap in the hits obtained"
                                             " from GWAS. It takes a full summary statistics file for a single trait"
                                             " and another file containing the GWAS hits from multiple traits. It then"
                                             " swaps in the hits from the GWAS hits file into the full summary"
                                             "statistics file. The resulting file can be used to generate a Manhattan"
                                             " plot that will show all of the hits across several traits.")
parser.add_argument("template_file",
                    metavar="TEMPLATE_FILEPATH",
                    help="Path to the template sumstats file.")
parser.add_argument("hits_file",
                    metavar="HITS_FILEPATH",
                    help="Path to the directory containing the GWAS hits files.")
parser.add_argument('-a', '--alias_file', metavar="ALIAS_FILEPATH",
                    help="Path to the file containing the aliases for the traits.")
parser.add_argument("-o", "--output_filepath",
                    help="Path to the directory where the output file will be written to.")

args = parser.parse_args()


hits_df = pl.read_csv(args.hits_file, sep='\t',
                      columns=["ID", "beta", "chi2", "pval", "Marker", "chromosome", "position", "OA", "EA",
                               "EAF", "Info", "phenotype"])

# If there is an alias file, read it and add the alias column to the hits_df
if args.alias_file:
    if not os.path.isfile(args.alias_file):
        raise ValueError('The provided alias file does not exist.')

    # Dictionary used to assign alias to each phenotype
    alias_dict = {}
    with open(args.alias_file, 'r') as in_file:
        for line in in_file:
            split_line = line.split(',')
            if len(split_line) == 2:
                alias_dict[split_line[0].replace(' ', '')] = split_line[1].rstrip()
            elif len(split_line) < 2:
                alias_dict[split_line[0].replace(' ', '')] = "Other"
            else:
                raise ValueError('The alias file is not formatted correctly. Each line must have 2 entries max.')

    # add column 'alias' to the hits_df dataframe using the alias_dict to map the values
    hits_df = hits_df.with_columns((pl.col('phenotype').apply(lambda x: alias_dict[x])).alias('alias'))

else:
    hits_df = hits_df.with_columns((pl.col('phenotype')).alias('alias'))  # If no alias file, just repeat the phenotype

template_df = pl.read_csv(args.template_file, sep="\t",
                          columns=["ID", "beta", "chi2", "pval", "Marker", "chromosome", "position", "OA", "EA",
                                   "EAF", "Info", "phenotype"])
template_df = template_df.with_columns((pl.col('phenotype')).alias('alias'))  # Add the alias column to the template_df

# Replace the relevant entries in the templateGWAS_df dataframe with the contents of the hits_df dataframe
out_df = pl.concat([template_df, hits_df])
# TODO: Ensure that the merge is working correctly. This is the most important step.

# Sort the created sumstats by chromosome and position before writing to file
# If chromosome names start with 'chr', remove it before sorting
if str(out_df.head(1).select('chromosome')[0, 0]).startswith('chr'):
    out_df.chromosome = out_df.with_column(pl.col('chromosome').str.replace('chr', ''))
# Drop duplicates (ignores the alias column)
out_df = out_df.unique(subset=['ID', 'beta', 'chi2', 'pval', 'Marker', 'chromosome', 'position', 'OA', 'EA', 'EAF',
                               'Info', 'phenotype'])
out_df = out_df.sort([pl.col('chromosome'), pl.col('position')])  # sort by column chromosome

# Write the output_manhattan_file
out_df.write_csv(args.output_filepath + '/combined_manhattan.txt', sep='\t')
