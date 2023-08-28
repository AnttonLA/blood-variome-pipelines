import os
import sys
import argparse
import re
import pandas as pd
import hail as hl

# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_file", required=True, help="Path to sample file")
parser.add_argument("-b", "--bgens_dir", required=True, help="Path to the directory containing the bgen files. BGEN"
                                                             "files need to be sorted and 8-")
parser.add_argument("-p", "--pheno_file", required=True, help="Path to phenotype file. It is recommended that the"
                                                              "phenotype values be Rank Inverse Normalized")
parser.add_argument("-x", "--sex_cov", help="File containing sex information for each sample. Used as a covariate")
parser.add_argument("-a", "--ancestry_cov", help="File containing PCA information for each sample. Used as a "
                                                 "covariate in place of sample donor ancestry")
parser.add_argument("-o", "--output_dir", required=True, help="Output folder where the results will be saved")

args = parser.parse_args()

# Check input arguments
for file in [args.sample_file, args.pheno_file]:
    if not os.path.isfile(file):
        raise ValueError(f"File {file} does not exist")

if not os.path.isdir(args.bgens_dir) or len(os.listdir(args.bgens_dir)) == 0:
    raise ValueError(f"BGEN folder {args.bgens_dir} does not exist or is empty")

if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)
    sys.stdout.write(f"Created output directory {args.output_dir}\n")


########################################################################################################################
# Utility functions
########################################################################################################################

def check_index(bgens_list: list) -> None:
    """
    Check if the BGEN files are indexed. If they are not, ask the user if they want to index them.

    :param bgens_list: List of absolute paths of the BGEN files
    """
    lacking_index = False
    for bgen in bgens_list:
        if not os.path.isdir(bgen + '.idx2'):
            lacking_index = True
            break

    # If any of the bgen files are not indexed, ask the user if they want to index them
    if lacking_index:
        print("At least some of the bgen files are not indexed. Do you want to index them? (y/n)")
        while True:
            answer = input()
            if answer == 'y':
                hl.index_bgen(bgens_list)
                break
            elif answer == 'n':
                print("Exiting")
                sys.exit(0)
            else:
                print("Please answer 'y' or 'n'")
    else:
        print("All bgen files are indexed. Proceeding...")


def sanitize_filename(filename: str) -> str:
    """
    Remove invalid characters from a filename. Used when we want to use a column from a dataframe as a filename to make
    sure we will get a valid name.
    """
    invalid_chars = r' <>:"/\|?*'
    valid_chars = [c for c in filename if c not in invalid_chars]
    sanitized_filename = ''.join(valid_chars).strip(' .')
    return sanitized_filename


def extract_chromosome_number(file_path: str):
    """Extract the chromosome number from a file path."""
    # Use regular expression to find the number after 'chr' and before '.bgen'
    match = re.search(r'chr(\d+)\.bgen', file_path)
    if match:
        return int(match.group(1))
    else:
        return None


def sort_bgen_files_by_chromosome(bgen_paths: list) -> list:
    """Sort a list of BGEN file paths by chromosome number."""
    return sorted(bgen_paths, key=extract_chromosome_number)


########################################################################################################################
# File loading and checking
########################################################################################################################

# Load the sample file. We don't need to do anything with it, but we can make sure it has the right dimensions.
# See format: https://www.cog-genomics.org/plink/2.0/formats#sample
sample_df = pd.read_csv(args.sample_file, sep=" ", header=0)

for col in ["ID_1", "ID_2", "missing", "sex"]:
    if col not in sample_df.columns.to_list():
        raise ValueError(f"Could not find column {col} in sample file {args.sample_file}. Please check the format "
                         f"of the sample file is correct: https://www.cog-genomics.org/plink/2.0/formats#sample")

print(f"Sample file loaded successfully. It contains {sample_df.shape[0]} samples.")

# Load the bgen files
bgens_abs_paths_list = []
for f in os.listdir(args.bgens_dir):
    if os.path.isfile(os.path.join(args.bgens_dir, f)) and f.endswith('.bgen'):
        bgens_abs_paths_list.append(os.path.join(args.bgens_dir, f))

# Sort the bgen files by chromosome number
bgens_abs_paths_list = sort_bgen_files_by_chromosome(bgens_abs_paths_list)

# Check if files are indexed. If not, offer option to index them
check_index(bgens_abs_paths_list)

# Check the phenotype file. It MUST contain a column named "Sample_ID" and at least one other column with a phenotype
pheno_df = pd.read_csv(args.pheno_file, sep="\t", header=0)
if "Sample_ID" not in pheno_df.columns.to_list():
    raise ValueError(f"Could not find column Sample_ID in phenotype file {args.pheno_file}. Please make sure a "
                     f"Sample_ID column is present in the phenotype file.")
if len(pheno_df.columns.to_list()) < 2:
    raise ValueError(f"Phenotype file {args.pheno_file} does not contain any phenotypes.")

phenotypes_list = pheno_df.columns.to_list()
phenotypes_list.remove("Sample_ID")
print(f"Phenotype file loaded successfully. It contains {pheno_df.shape[0]} samples and {len(phenotypes_list)} "
      f"phenotype(s).")

# Check the sample IDs shared by the sample file and the phenotype file
shared_sample_ids = set(sample_df["ID_1"].to_list())
shared_sample_ids.intersection_update(set(pheno_df["Sample_ID"].to_list()))
if len(shared_sample_ids) == 0:
    raise ValueError(f"No sample IDs are shared between the sample file and the phenotype file. Please make sure "
                     f"both files contain the same sample IDs.")
else:
    print(f"{len(shared_sample_ids)} sample IDs are shared between the sample file and the phenotype file.")

if args.sex_cov is not None:
    sex_df = pd.read_csv(args.sex_cov, sep="\t", header=0)
    # Check that the df has two columns only
    if len(sex_df.columns.to_list()) != 2:
        raise ValueError("Sex covariate information file must be a tab separated file with two columns. Provided file "
                         f"'{args.sex_cov}' was not valid. Please provide a valid file.")
    if 'Sample_ID' not in sex_df.columns.to_list():
        raise ValueError("Sex covariate file must contain a column named 'Sample_ID'. Please provide a valid file.")

if args.ancestry_cov is not None:
    ancestry_df = pd.read_csv(args.ancestry_cov, sep="\t", header=0)
    # Check that the df has at least two columns
    if len(ancestry_df.columns.to_list()) < 2:
        raise ValueError("Ancestry covariate information file must be a tab separated file with at least two columns. "
                         f" Provided file '{args.ancestry_cov}' was not valid. Please provide a valid file.")
    if 'Sample_ID' not in ancestry_df.columns.to_list():
        raise ValueError("Ancestry covariate file must contain a column named 'Sample_ID'. Please provide a valid file.")

print(sex_df)
sys.exit()
########################################################################################################################
# GWAS
########################################################################################################################

hl.init()

# Load the bgen files. Generate Matrix Table from BGENs + sample file
mt = hl.import_bgen(bgens_abs_paths_list, entry_fields=['GT', 'GP'], sample_file=args.sample_file)

# Add row parameter 'variant_qc' to MatrixTable
mt = hl.variant_qc(mt)

# Minimal QC and data filtering
# Filter out variants with MAF < 0.05
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.05)
# HW equilibrium filtering
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)

# Annotate the matrix table with phenotype information
table = hl.import_table(args.pheno_file, delimiter='\t', missing="", impute=True,
                        types={'Sample_ID': hl.tstr}).key_by('Sample_ID')
mt = mt.annotate_cols(**table[mt.s])

# Run GWAS for each phenotype
for i, phenotype in enumerate(phenotypes_list):
    print(f"Commencing GWAS {i + 1}/{len(phenotypes_list)}.\n\tPhenotype: {phenotype}\n")
    gwas = hl.linear_regression_rows(y=mt[phenotype], x=mt.GT.n_alt_alleles(), covariates=[1.0],
                                     pass_through=[mt.rsid, mt.variant_qc])
    gwas = gwas.key_by('rsid')
    results = gwas.select(chromosome=gwas.locus.contig, position=gwas.locus.position, OA=gwas.alleles[0],
                          EA=gwas.alleles[1], EAF=gwas.variant_qc.AF[1], pval=gwas.p_value, beta=gwas.beta,
                          tstat=gwas.t_stat, n=gwas.n)

    # Save results to file
    output_file = os.path.join(args.output_dir, f"GWAS_{sanitize_filename(phenotype)}.tsv")
    results.export(output_file, header=True)

    inflation = hl.methods.lambda_gc(gwas.p_value)
    print(f'GWAS #{i + 1} complete!\n\tInflation factor: {inflation}\n')

    if i > 2:
        break  # TODO: Remove this line when testing is complete
