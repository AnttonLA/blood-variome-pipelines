# Causal Gene Assignment Pipeline

This pipeline is designed to assign putative causal genes to a list of genetic variants.
It uses data from [ImmuNexUT](https://www.immunexut.org/top) and the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)
to find genes that have their expression affected by the variants.

The pipeline makes heavy use of [GORpipe](https://docs.gorpipe.org/index.html), a tool for processing genomic data.
It is written in [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system.

## Running the pipeline

To run this pipeline, you will need GORpipe and Snakemake installed, as well as data from the eQTL Catalogue and ImmuNexUT
that has been processed by GORpipe.

Start by editing the `config.yaml` file to point to the correct directories, data files and executables.
Then, run the pipeline with:

    snakemake --cores all
