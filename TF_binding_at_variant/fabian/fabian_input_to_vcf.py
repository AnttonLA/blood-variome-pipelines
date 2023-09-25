import os
import argparse

"""
This script takes a list of variants in the format "chr:posOA>EA" and converts it to a VCF file.
This is necessary to run FABIAN-Variant for more than 100 variants.
The maximum amount of variants that can be run at once is 10k. If you want to run more than that, you need to split the
list into smaller chunks and run them separately.
"""

parser = argparse.ArgumentParser(description="Convert a variant list to the VCF format required by FABIAN-Variant")
parser.add_argument("input_file", help="Path to input variant file")
parser.add_argument("-o", "--output_dir", required=True, help="Path to output directory")

args = parser.parse_args()

# Make sure the output directory exists
if not os.path.isdir(args.output_dir):
    raise ValueError(f"Output directory {args.output_dir} does not exist")


# List of variants in Fabian input format
with open(args.input_file, "r") as input_file:
    variants = input_file.readlines()

num_variants = len(variants)
if num_variants > 10000:
    num_output_files = num_variants // 10000 + 1
    print(f"Warning: FABIAN-Variant can only run 10k variants at a time. You have {num_variants} variants. "
          f"The output will be split into {num_output_files} files.")
else:
    num_output_files = 1


# Create a VCF file for each 10k variant chunk
for i in range(num_output_files):
    vcf_file_name = f"FABIAN_INPUT_{i+1}.vcf"
    with open(os.path.join(args.output_dir, vcf_file_name), "w") as vcf_file:
        # Define the header for the VCF file
        header = """##fileformat=VCFv4.2
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001"""

        vcf_file.write(header)

        for variant in variants[i*10000:(i+1)*10000]:
            parts = variant.split(":")
            chrom = parts[0]
            # Split parts[1] into pos, ref, alt. First take all the numbers from the beginning of the string
            pos = ""
            ref = ""
            alt = ""
            for char in parts[1]:
                if char.isdigit():
                    pos += char
                else:
                    break
            # Then take the rest of the string
            ref_alt = parts[1].replace(pos, "")
            ref = ref_alt.split(">")[0]
            alt = ref_alt.split(">")[1].replace("\n", "")

            vcf_line = f"{chrom.lstrip('chr')}\t{pos}\t{variant.rstrip()}\t{ref}\t{alt}\t2000\t.\t.\tGT:DP\t0/1:154\n"
            vcf_file.write(vcf_line)

    print(f"VCF file {vcf_file_name} created at {args.output_dir}")
