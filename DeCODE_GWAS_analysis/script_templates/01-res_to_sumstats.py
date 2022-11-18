import pandas as pd
from scipy.stats import chi2
import os

"""This script takes a DeCODE .rse file (table with ID, beta and chi-square columns) and filters the table by pval threshold.
"""

files_filepath = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio"
variant_info_file = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/variant_info.txt"

output_filepath = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/combined_output/all_significant_variant_tables"
#gwas_file = "SWE_Swedes_Blood_variome_Bpanel_CD16negCD56posNK_div_Bpanel_CD45pos_Frequency_adjSexPhaCohPC_InvNorm_12102022.res"

# Read variant infor table
f2 = pd.read_csv(variant_info_file, sep='\t', chunksize=10 ** 6)
list_of_f2_chunks = []
for chunk in f2:
    list_of_f2_chunks.append(chunk)
variant_info_df = pd.concat(list_of_f2_chunks)
print("Variant info table loaded!\nProceeding with result filtering...")

for i, gwas_file in enumerate(os.listdir(files_filepath)):
    if gwas_file.endswith('.res'):
        print(f"File {i + 1} of {len(os.listdir(files_filepath))}")
        #f = pd.read_csv(files_filepath + '/' + gwas_file, sep=" ", chunksize=10**6)
        f = pd.read_csv(files_filepath + '/' + gwas_file, delim_whitespace=True, chunksize=10**6)

        list_of_downsampled_chunks = []
        print(f"Reading file {gwas_file}")
        for chunk in f:
            chunk.columns = ['ID', 'beta', 'chi2']

            # Assuming DF=1, chi2 > 23 is roughly pval < 1E-6
            list_of_downsampled_chunks.append(chunk[chunk.chi2 > 23].copy())  # Take

        phenotype = '_'.join(gwas_file.rstrip('.txt').split('_')[4:-3])  # Take pheno name from file name
        concat_df = pd.concat(list_of_downsampled_chunks)
        concat_df = concat_df.assign(phenotype=phenotype)
        concat_df.insert(3, 'pval', concat_df.chi2.apply(lambda x: chi2.sf(x, 1)))  # Convert chi2 to pval

        # Merge the two dfs by ID
        merged_df = concat_df.merge(variant_info_df, on='ID', how='left')
        merged_df = merged_df[["ID", "beta", "chi2", "pval", "Marker", "OA", "EA", "EAF", "Info", "phenotype"]].copy()

        merged_df.insert(5, 'chromosome', merged_df.Marker.apply(lambda x: x.split(':')[0]))  # Extract chr from Marker
        merged_df.insert(6, 'position', merged_df.Marker.apply(lambda x: x.split(':')[1]))  # Extract pos from Marker

        merged_df.to_csv(output_filepath + '/' + phenotype + '_significant_variants.txt', sep='\t', index=False)
