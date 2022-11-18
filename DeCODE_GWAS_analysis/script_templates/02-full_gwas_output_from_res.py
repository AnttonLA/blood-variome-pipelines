import pandas as pd
from scipy.stats import chi2
import os

"""This script generates a full "template Manhattan" file from the DeCODE data to then replace the significant variants into"""

output_filepath = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/combined_output/"

files_filepath = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio"
variant_info_file = "/home/antton/Projects/Immune_GWAS/data/raw/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/variant_info_extended.txt"

gwas_file = "SWE_Swedes_Blood_variome_Bpanel_CD16negCD56posNK_div_Bpanel_CD45pos_Frequency_adjSexPhaCohPC_InvNorm_12102022.res"

# Read variant infor table
f2 = pd.read_csv(variant_info_file, sep='\t', chunksize=10 ** 6)
list_of_f2_chunks = []
for chunk in f2:
    list_of_f2_chunks.append(chunk)
variant_info_df = pd.concat(list_of_f2_chunks)
print("Variant info table loaded!\nProceeding with chi2 score file")

if gwas_file.endswith('.res'):
    f = pd.read_csv(files_filepath + '/' + gwas_file, delim_whitespace=True, chunksize=10**6)

    list_of_downsampled_chunks = []
    print(f"Reading file {gwas_file}")
    for i, chunk in enumerate(f):
        print(f"Chunk {i} loaded")
        chunk.columns = ['ID', 'beta', 'chi2']

        # Assuming DF=1, chi2 > 23 is roughly pval < 1E-6
        list_of_downsampled_chunks.append(chunk[chunk.chi2 > 2].copy())  # append full chunk to list

    phenotype = '_'.join(gwas_file.rstrip('.txt').split('_')[4:-3])  # Take pheno name from file name
    concat_df = pd.concat(list_of_downsampled_chunks)
    concat_df = concat_df.assign(phenotype=phenotype)
    print("Conveting chi2 to pval")
    concat_df.insert(3, 'pval', concat_df.chi2.apply(lambda x: chi2.sf(x, 1)))  # Convert chi2 to pval

    print("Commencing merge with variant info table")
    # Merge the two dfs by ID
    merged_df = concat_df.merge(variant_info_df, on='ID', how='left')
    merged_df = merged_df[["ID", "beta", "chi2", "pval", "Marker", "Chr", "Pos", "OA", "EA", "EAF", "Info", "phenotype"]].copy()

    merged_df.rename(columns={"Chr": "chromosome", "Pos": "position"}, inplace=True)
    print("Done! Saving to file:")
    merged_df.to_csv(output_filepath + '/template_manhattan.txt', sep='\t', index=False)
