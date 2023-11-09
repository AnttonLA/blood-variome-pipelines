import os
import sys
import polars as pl
import argparse

"""
Convert a variant list into the format required by FABIAN-Variant.
"""


def variant_list_to_fabian_input_vcf(variant_file, output_dir):
    """
    This function will take a file containing information about variants and convert it to a VCF file.
    The resulting file will be used as an input for FABIAN-Variant.
    
    :param variant_file: 
    :param output_dir: 
    :return: 
    """
    # Check if the input file exists
    if not os.path.isfile(variant_file) or os.path.getsize(variant_file) == 0:
        raise ValueError(f"Input variant list file {variant_file} does not exist or is empty")

    # Ensure the input file has the necessary columns (ID, Chrom, Pos, OA, EA)
    df = pl.read_csv(variant_file, separator='\t', has_header=True,
                     dtypes={"ID": pl.Utf8, "Chrom": pl.Utf8, "Pos": pl.Int64, "OA": pl.Utf8, "EA": pl.Utf8})
    if not {"ID", "Chrom", "Pos", "OA", "EA"}.issubset(set(df.columns)):
        raise ValueError(
            f"Input variant list file {variant_file} does not have the necessary columns. Please ensure "
            f"the file contains the columns ID, Chrom, Pos, OA, EA.")

    # Remove the 'chr' prefix from the Chrom column
    df = df.with_columns(pl.col("Chrom").str.replace("chr", "").alias("Chrom"))
    # Change its dtype to int
    df = df.with_columns(pl.col("Chrom").cast(pl.Int64))

    # Rename the columns to match the VCF format
    df = df.rename({"Chrom": "#CHROM", "Pos": "POS", "OA": "REF", "EA": "ALT"})

    # Add the missing columns
    df = df.with_columns([pl.lit(100).alias("QUAL"), pl.lit(".").alias("FILTER"), pl.lit(".").alias("INFO"),
                          pl.lit("GT:DP").alias("FORMAT"), pl.lit("0/1:154").alias("NA00001")])

    # Reorder the columns and sort by chrom and pos
    df = df.select(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA00001"])
    df = df.sort(["#CHROM", "POS"])

    # If df is too large the output will need to be split into several VCF files
    num_variants = df.shape[0]
    if num_variants > 10000:
        num_output_files = num_variants // 10000 + 1
        sys.stdout.write(
            f"\nWARNING: FABIAN-Variant can only run 10k variants at a time. You have {num_variants} variants."
            f" The output will be split into {num_output_files} files.")
    else:
        num_output_files = 1

    # Create a VCF file for each 10k variant chunk
    for idx, sub_df in enumerate(df.iter_slices()):
        vcf_file_name = f"FABIAN_INPUT_{idx + 1}.vcf"
        sub_df.write_csv(os.path.join(output_dir, vcf_file_name), separator="\t", has_header=True)
        print(f"VCF file {vcf_file_name} created at {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a variant list to the format required by FABIAN-Variant")
    parser.add_argument("input_file", help="Path to input variant list file")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Path to output directory where the VCF file(s) will be "
                             "created")

    args = parser.parse_args()

    variant_list_to_fabian_input_vcf(args.input_file, args.output_dir)
