#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script that runs the perfectos-ape java ape.jar for a given snplist
The snplist file contains a list of SNPs in the following format:
ID CHR POS OA EA

Original author: Ludvig Ekdahl
Modified by: Antton Lamarca
"""
import sys
import os
import argparse


def create_hocomocco_input_file(snplist: str, hg38fa: str, samtools: str, verbose: bool = False) -> None:
    """Parse input file and extract chromosome, position. Then expand every position to a total of 60 bp with the SNP
    in the middle. The SNP itself is written in brackets (e.g. [A/T]). The output is written to a tmp file and later
    used as input for the ape.jar. Indels are ignored and saved to a separate list that is outputted at the end.

    :param snplist: Path to the snplist file
    :param hg38fa: Path to the hg38 fasta file
    :param samtools: Path to the samtools executable
    :param verbose: If True, print more information to stdout/stderr
    :return:
    """
    # Input checks
    if not os.path.isfile(snplist):
        raise ValueError("SNP list file does not exist")
    if not os.path.isfile(hg38fa):
        raise ValueError("hg38 fasta file does not exist")

    # Iterate through snplist file and create temporary file with input for ape.jar for each SNP.
    # Indels are stored in another temp file, but nothing is done with them.
    indels = []
    with open(snplist, "r") as snp_file:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        with open(f"{script_dir}/tmp/tmp.txt", "w") as tmp_file:
            # Skip header
            next(snp_file)
            i = 0
            for line in snp_file:
                i += 1
                line = line.strip().split()
                if len(line) != 5:
                    raise ValueError(f"SNP list file is not in the correct format. Error in line {i}: {line}")
                snp_id, chrom, pos, oa, ea = line[0], line[1], line[2], line[3], line[4]

                # TODO : Manage indels
                # EXAMPLE OF INDEL THAT CURRENTLY FAILS HERE
                # chr19:49525049_CCTGCTGCTGCTG_CCTGCTG	TCCTCCCTGAGTCTGACCATCTTCCAT[CCTGCTGCTGCTG/CCTGCTG]ctgctgctgctgctgctgcGG
                # THIS IS A DELETION OF 4 BASES
                # THE OA CONTAINS ENOUGH BASES TO MAP IT TO THE REFERENCE
                # THE BRACKET SECTION SHOULD REPLACE THE ENTIRE SECTION OF THE REFERENCE THAT 'OA' MAPS TO
                # I.E collect ref -> MATCH CCTGCTGCTGCTG -> REPLACE WITH [CCTGCTGCTGCTG/CCTGCTG] AND THEN CONTINUE SEQUENCE
                # - Ludvig
                # TODO: FOR NOW SKIP, LATER IMPLEMENT THESE AS DIFFERENT FUNCTION

                # Collect indels for later processing, skip for now
                if len(oa) > 1 or len(ea) > 1:
                    if verbose:
                        sys.stdout.write(f"Skipping indel {snp_id}\n")
                    indels.append(snp_id)
                    continue

                flank = int((30 - len(ea) / 2) + 0.5)
                start = int(pos) - flank
                end = int(pos) + flank
                # Now pull out the sequence from the fasta file
                # The fasta file is indexed with samtools faidx
                samtoolsfaidx_return = os.popen(f"{samtools} faidx {hg38fa} {chrom}:{start}-{end}").read().split("\n")

                # The sequence can vary in length between lines. Merge all lines beyond the first. Remove "\n".
                actual_sequence = "" + "".join(samtoolsfaidx_return[1:]).strip().replace("\n", "")

                # We replace the SNP with OA/EA inside hard brackets in the sequence at start + flank + 1
                snp_sequence = actual_sequence[:flank] + "[" + oa + "/" + ea + "]" + actual_sequence[flank + len(ea):]

                # Construct a perfectos-APE input line, which is snp_ID "\t" snp_sequence. Write to tmp file.
                perfectos_ape_input_line = f"{snp_id} {snp_sequence}\n"
                tmp_file.write(perfectos_ape_input_line)

        # Write indels to a separate file
        with open(f"{script_dir}/tmp/indels.txt", "w") as indel_file:
            for indel in indels:
                indel_file.write(indel + "\n")


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--snplist", required=True,
                        help="Path to snplist file, see snps_for_motif_ZA.txt for example")
    parser.add_argument("-r", "--ref", required=True,
                        help="Path to reference fasta file. Needs to have been indexed with samtools faidx")
    parser.add_argument("--samtools", required=False, default="samtools",
                        help="Path to samtools binary. Not needed if samtools is in PATH (default: %(default)s)")
    parser.add_argument("--hocomoco", required=False, default="./pwm/hocomoco_11_human",
                        help="Path to hocomoco database (default: %(default)s)")
    parser.add_argument("-v", "--verbose", required=False, action="store_true",
                        help="Print more information to stdout/stderr")
    parser.add_argument("-o", "--output", required=False, default="./results/results.txt",
                        help="Path to output file. (default: %(default)s)")

    args = vars(parser.parse_args())

    script_dir = os.path.dirname(os.path.realpath(__file__))  # Directory where this script is located
    ape_jar = f"{script_dir}/ape.jar"  # Path to the perfectos-ape jar file

    # Check that the .jar file exists
    if not os.path.isfile(ape_jar):
        raise FileNotFoundError(f"Could not find ape.jar in expected location {script_dir}")

    snplist = args["snplist"]
    hocomoco_db = args["hocomoco"]
    hg38fa = args["ref"]  # Path to the hg38 reference genome fasta file
    samtools = args["samtools"]

    # Pull out ref sequences from reference fasta file and create a PERFECTOS-APE input file
    create_hocomocco_input_file(snplist, hg38fa, samtools)

    sys.stdout.write(f"Created perfectos-ape input file, running PERFECTOS-APE... (this may take a while)\n"
                     f"Output will be continuously written to {args['output']}\n")

    # Call the perfectos-ape java ape.jar with the given snplist
    # Catch standard output and redirect it to a results file
    os.system(f"java -jar {ape_jar} {hocomoco_db} {script_dir}/tmp/tmp.txt > {args['output']}")

    sys.stdout.write("Done!\n")
