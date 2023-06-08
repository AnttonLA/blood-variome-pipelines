#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script that runs the perfectos-ape java ape.jar for a given snplist
The snplist file contains a list of SNPs in the following format:
ID CHR POS OA EA

Original author: Ludvig Ekdahl
Modified by: Antton Lamarca
"""
import os
import argparse


def create_hocomocco_input_file(snplist, hg38fa, samtools) -> None:
    """Parse input file and extract chromosome, position. Then expand every position to a total of 60 bp with the SNP
    in the middle. The SNP itself is written in hard brackets like [A/T]. The output is written to a new tmp file and
    then used as input for the ape.jar. Indels are ignored and saved to a separate list that is outputted at the end.

    :param snplist: Path to the snplist file
    :param hg38fa: Path to the hg38 fasta file
    :param samtools: Path to the samtools executable
    :return:
    """

    indels = []
    with open(snplist, "r") as snp_file:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        with open(f"{script_dir}/tmp/tmp.txt", "w") as tmp_file:
            # Skip header
            next(snp_file)
            n = 0
            for line in snp_file:
                line = line.strip().split()
                snp_id = line[0]
                chrom = line[1]
                pos = line[2]
                oa = line[3]
                ea = line[4]
                # check length of ea allele, adjust length of flanking sequence
                # TODO : SHOULD TOOL EVEN ALLOW INDEL
                # EXAMPLE OF INDEL THAT CURRENTLY FAILS HERE
                # chr19:49525049_CCTGCTGCTGCTG_CCTGCTG	TCCTCCCTGAGTCTGACCATCTTCCAT[CCTGCTGCTGCTG/CCTGCTG]ctgctgctgctgctgctgcGG
                # THIS IS A DELETION OF 4 BASES
                # THE OA CONTAINS ENOUGH BASES TO MAP IT TO THE REFERENCE
                # THE BRACKET SECTION SHOULD REPLACE THE ENTIRE SECTION OF THE REFERENCE THAT 'OA' MAPS TO
                # I.E collect ref -> MATCH CCTGCTGCTGCTG -> REPLACE WITH [CCTGCTGCTGCTG/CCTGCTG] AND THEN CONTINUE SEQUENCE
                # TODO: FOR NOW SKIP, LATER IMPLEMENT THESE AS DIFFERENT FUNCTION
                if len(oa) > 1 or len(ea) > 1:
                    # report indel
                    print(f"Skipping indel {snp_id}")
                    indels.append(snp_id)
                    continue
                flank = int((30 - len(ea) / 2) + 0.5)
                start = int(pos) - flank
                end = int(pos) + flank
                # Now pull out the sequence from the fasta file
                # The fasta file is indexed with samtools faidx
                samtoolsfaidx_return = os.popen(f"{samtools} faidx {hg38fa} {chrom}:{start}-{end}").read().split("\n")

                # The fasta sequence can be variable length per line, so we merge all lines beyond the first
                # And remove the newline characters
                actual_sequence = "" + "".join(samtoolsfaidx_return[1:]).strip().replace("\n", "")

                # We replace the SNP with OA/EA inside hard brackets in the sequence at start + flank + 1
                snp_sequence = actual_sequence[:flank] + "[" + oa + "/" + ea + "]" + actual_sequence[flank + len(ea):]

                # Construct a perfectos-APE input line, which is snp_ID "\t" snp_sequence
                perfectos_ape_input_line = snp_id + " " + snp_sequence + "\n"

                # Write the perfectos-ape input line to the tmp file
                tmp_file.write(perfectos_ape_input_line)

        # Write indels to a separate file
        with open(f"{script_dir}/tmp/indels.txt", "w") as indel_file:
            for indel in indels:
                indel_file.write(indel + "\n")


if __name__ == "__main__":
    # Parse arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--snplist", required=True,
                    help="(REQUIRED) Path to snplist file, see snps_for_motif_ZA.txt for example")
    ap.add_argument("--ref", required=True, default="./hg38.fa",
                    help="(REQUIRED) Path to reference fasta file, indexed with samtools faidx")
    ap.add_argument("--samtools", required=False, default="samtools",
                    help="Path to samtools binary, not required if samtools is in path")
    ap.add_argument("--hocomoco", required=False, default="./pwm/hocomoco_11_human",
                    help="Path to hocomoco database")
    ap.add_argument("-o", "--output", required=False, default="./results/results.txt",
                    help="Path to output file, default is ./results/results.txt")

    args = vars(ap.parse_args())

    script_dir = os.path.dirname(os.path.realpath(__file__))  # Directory where this script is located
    # Path to the perfectos-ape jar file
    ape_jar = f"{script_dir}/ape.jar"

    # Path to the snplist file
    snplist = args["snplist"]

    # Path to the hocomoco database
    hocomoco_db = args["hocomoco"]

    # Path to the hg38 fasta file
    hg38fa = args["ref"]  # "/media/ludvig/cbio3/projects/Zain_2021/hg38FASTA/hg38.fa"

    # path to samtools binary
    samtools = args["samtools"]  # "/home/ludvig/Programs/samtools/samtools/samtools"

    # This pulls out ref sequences from reference fasta file and creates a perfectos-ape input file
    create_hocomocco_input_file(snplist, hg38fa, samtools)

    # report that the perfectos-ape input file has been created
    print("Created perfectos-ape input file, running perfectos-ape... (this may take a while)\nOutput will be "
          f"continuously written to {args['output']}")

    # Call the perfectos-ape java ape.jar with the given snplist
    # Catch stanard output and redirect it to a results file
    os.system(f"java -jar {ape_jar} {hocomoco_db} {script_dir}/tmp/tmp.txt > {args['output']}")

    # report that the perfectos-ape run has been completed
    print("Completed")
