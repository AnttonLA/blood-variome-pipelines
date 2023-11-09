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

## References

Fayrouz Hammal, Pierre De Langen, Aurélie Bergon, Fabrice Lopez, Benoit Ballester, **ReMap 2022: a database of Human, Mouse, 
Drosophila and Arabidopsis regulatory regions from an integrative analysis of DNA-binding sequencing experiments**, 
*Nucleic Acids Research*, Volume 50, Issue D1, 7 January 2022 

Robin Steinhaus, Peter N Robinson, Dominik Seelow, **FABIAN-variant: predicting the effects of DNA variants on
transcription factor binding**, *Nucleic Acids Research*, Volume 50, Issue W1, 5 July 2022, Pages W322–W329,
https://doi.org/10.1093/nar/gkac393
