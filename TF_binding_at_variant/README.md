# Transcription Factor Assignment Pipeline

This is a pipeline intended to find the Transcription Factors (TFs) that are most likely being affected by a given
non-coding genetic variant.

It takes as input a file with the following format:

```
ID	Chrom	Pos	OA	EA
variant_id	chr1	23935190	G	T
```

You will also need:
* ReMap metadata BED file, bgzipped and tabix-indexed

## Steps
Currently, there are only three steps to this pipeline:

1. Look up entries of the ReMAP database (https://remap2022.univ-amu.fr/) to find Transcription Factors that
bind over a given variant.
2. Use either FABIAN-Variant (https://www.genecascade.org/fabian/) to find the transcription factors whose binding motif is likely
 disrupted (or created/strengthened) by the variant.
3. Finally, the pipeline will then assign those TFs that appear in **BOTH** to the variant.

It **requires tabix** to perform the lookup on the ReMap metadata table and find all the relevant studies. 
It also **requires samtools** to index the BAM files from the ReMap database.
You can find a more detailed description of the pipeline below.

## Installation

The pipeline is written in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses
[conda](https://docs.conda.io/en/latest/) to manage the environments.

It is not available as a standalone pipeline, so you will need to download the entire repository.  
You can do so with:

```
git clone https://github.com/AnttonLA/blood-variome-pipelines.git
```

Then, create the conda environment with:

```
cd blood-variome-pipelines/TF_binding_at_variant
conda env create -f environment.yaml
```

## Usage

Edit the paths in the config.yaml file so that they point to the correct locations in your system. 
Then, run the pipeline by running Snakemake with:

```
snakemake --cores=all
```

### Future Additions

Some steps that I would like to add to the pipeline include:

- [ ] Use [Selenium](https://www.selenium.dev/) to automate uploading the files to FABIAN-Variant.
- [ ] Incorporate [ffq](https://github.com/pachterlab/ffq) to download the relevant data from ReMap. 


## Description of the pipeline

### Step 1: Find TFs that bind over the variants
The pipeline will first look up the positions of the variants in the ReMap metadata table.
It will then find all the studies that contain a TF that binds over the variants and output a list with all the matches.
Finally, it will produce a second file where only the biotypes of interest are kept.

### Step 2: Find TFs whose binding motif is likely disrupted by the variant
The pipeline will then use FABIAN-Variant to find the TFs whose binding motif is likely disrupted by the variant.
The pipeline will create a VCF file that can be used as an input for FABIAN-Variant.

    NOTE: Currently the use of FABIAN-Variant is not automated.
    You will need to upload the VCF file to the FABIAN-Variant and download the results manually.

The pipeline will then parse the results from FABIAN-Variant. The full output file ("the data file") will be kept for 
future reference. The summary file ("the table file") is the one that will be used for the rest of the analysis, since
it is the one that contains the scores that have been averaged across models.

### Step 3: Assign TFs to the variant
Finally, the pipeline will assign the TFs that appear in **BOTH** the ReMap and FABIAN-Variant results to each variant.
Keep in mind that these results will vary considerably depending on the selected Biotypes and FABIAN S-score threshold.
## References

Fayrouz Hammal, Pierre De Langen, Aurélie Bergon, Fabrice Lopez, Benoit Ballester, **ReMap 2022: a database of Human, Mouse, 
Drosophila and Arabidopsis regulatory regions from an integrative analysis of DNA-binding sequencing experiments**, 
*Nucleic Acids Research*, Volume 50, Issue D1, 7 January 2022 

Robin Steinhaus, Peter N Robinson, Dominik Seelow, **FABIAN-variant: predicting the effects of DNA variants on
transcription factor binding**, *Nucleic Acids Research*, Volume 50, Issue W1, 5 July 2022, Pages W322–W329,
https://doi.org/10.1093/nar/gkac393
