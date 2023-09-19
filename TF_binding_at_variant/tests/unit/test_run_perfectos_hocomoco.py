import os
from perfectos_ape.run_perfectos_hocomoco import create_hocomoco_input_file


def test_create_hocomoco_input_file():
    snplist = "test_data/dummy_variant_list.tsv"
    reference_fasta = "/media/antton/cbio3/projects/Zain_2021/hg38FASTA/hg38.fa"

    create_hocomoco_input_file(snplist, reference_fasta, "samtools", "./tmp", True)

    # Assert that indels.txt and perfectos-ape_input.txt exist
    assert os.path.isfile("./tmp/indels.txt")
    assert os.path.isfile("./tmp/perfectos-ape_input.txt")

    # Remove both files
    os.remove("./tmp/indels.txt")
    os.remove("./tmp/perfectos-ape_input.txt")
