import pandas as pd
import sys
import os

"""
This script takes a list of rsIDs and checks the eQTL summary stats files for rows where these rsIDs appear.
Each eQTL data file is a different cell type in the ImmuNexUT dataset. Note that these files only contain associations
with FDR < 0.05. If we want to use other eQTL data files, we will need to implement some sort of p-value/FDR threshold.

Usage: python immunexut_eqtl_lookup.py <path_to_eQTL_data_folder> <list_of_rsIDs> <output_file_location> 

author: Antton Lamarca
2022-06-21
"""


def get_matching_rows(filepath, rsids, rsid_col_name='rsid', gene_col_name='gene', chr_col_name='chr',
                      pos_col_name='pos', pval_col_name='pval', beta_col_name='beta'):
    """
    This function takes a filepath to an eQTL summary stats file and a list of rsIDs. It returns a pandas dataframe with
    the eQTL summary stats for the requested rsIDs.

    :param filepath: global path to the eQTL summary stats file
    :param rsids: list of rsIDs to look for in the eQTL summary stats file
    :param rsid_col_name: name of the column in the eQTL summary stats file that contains the rsIDs
    :param gene_col_name: name of the column in the eQTL summary stats file that contains the gene names
    :param chr_col_name: name of the column in the eQTL summary stats file that contains the chromosome names
    :param pos_col_name: name of the column in the eQTL summary stats file that contains the position of the SNP
    :param pval_col_name: name of the column in the eQTL summary stats file that contains the p-value of the association
    :param beta_col_name: name of the column in the eQTL summary stats file that contains the beta value of the
    association
    :return: pandas DataFrame of rows where the rsIDs are present in the eQTL summary stats file
    """

    eQTL_df = pd.read_csv(filepath, sep='\t')  # load eQTL summary stats file into a df
    rsid_rows_df = eQTL_df[eQTL_df[rsid_col_name].isin(rsids)]  # get rows where rsIDs are in the eQTL data file
    # Take only columns we want to keep
    out_df = rsid_rows_df[[rsid_col_name, chr_col_name, pos_col_name,
                           gene_col_name, pval_col_name, beta_col_name]].copy()
    # Rename columns
    out_df.rename(columns={rsid_col_name: "rsid",
                           gene_col_name: "gene",
                           chr_col_name: "chromosome",
                           pos_col_name: "position_hg38",
                           pval_col_name: "pval",
                           beta_col_name: "slope"}, inplace=True)

    return out_df


if __name__ == '__main__':

    args = sys.argv[1:]

    if len(args) != 3:
        print("Wrong number of input arguments.")
        print("Usage: python immunexut_eqtl_lookup.py <path_to_eqtl_data_folder> <rsID_list> <output_file>")
        sys.exit(1)

    # Store path to eQTL data folder
    eqtl_folder = args[0]
    if not os.path.isdir(eqtl_folder):  # make sure 'eqtl_folder' is a path to a directory
        print("Path to eQTL data folder is not a directory.")
        sys.exit(1)
    if eqtl_folder.endswith('/'):
        eqtl_folder = eqtl_folder[:-1]

    # Read rsID list
    rsID_list_file = args[1]
    with open(rsID_list_file, 'r') as f:
        rsID_list = f.read().splitlines()

    # Save output file name
    output_filename = args[2]

    # Iterate over all eQTL data files
    df = pd.DataFrame()
    num_files = len([f for f in os.listdir(eqtl_folder) if f.endswith('.txt')])  # count # of .txt files in eqtl_folder

    counter = 0
    for file in os.listdir(eqtl_folder):
        if file.endswith('.txt'):
            counter += 1
            print(f'File {counter}/{num_files}: {file}')

            full_filepath = f'{eqtl_folder}/{file}'  # eQTL data global file path
            matching_rows_df = get_matching_rows(full_filepath, rsID_list, rsid_col_name="Variant_ID",
                                                 gene_col_name="Gene_name", chr_col_name="Variant_CHR",
                                                 pos_col_name="Variant_position_start",
                                                 pval_col_name="Forward_nominal_P",
                                                 beta_col_name="Forward_slope")
            matching_rows_df.insert(6, 'cell_type', '_'.join(file.split('_')[:-3]))  # Get cell type from file name

            df = pd.concat([df, matching_rows_df])

    # Sort by chromosome, position and gene name.
    # We temporarily add a chr_num column to the df to sort by chromosome number correctly, and drop de column after.
    df['chr_num'] = df['chromosome'].apply(lambda x: int(x.split('chr')[1]))
    df.sort_values(by=['chr_num', 'position_hg38', 'gene'], inplace=True)
    df.drop(columns=['chr_num'], inplace=True)
    df.to_csv(output_filename, sep='\t', index=False)
