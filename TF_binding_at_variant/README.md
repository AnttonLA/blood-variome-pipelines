Pipeline to find the Transcription Factor that is likely affected by a given variant.

It takes as input a file with the following format:

```
ID	Chrom	Pos	OA	EA
chr1:23935190_G_T	chr1	23935190	G	T
```

Currently, there are only two steps:

1. Look up entries of the ReMAP database (https://remap2022.univ-amu.fr/) to find Transcription Factors that
bind over a given variant.
2. Use PERFECTOS-APE to find the transcription factors whose binding motif likely disrupted by the variant.

The pipeline will then assign those TFs that appear in **BOTH** to the variant.

It **requires tabix** to perform the lookup on the ReMap metadata table and find all the relevant studies. 

### Usage

In order to run the pipeline, you will need:
* Input file with the format described above
* ReMap metadata BED file, bgzipped and tabix-indexed
* hg38 reference genome fasta file (used by PERFECTOS-APE)

You will also need to have samtools installed.

Edit the paths in the config.yaml file so that they point to the correct locations. 
Then, run the pipeline with:

```
snakemake --cores=all
```

### Future Additions

It would be nice to incorporate [ffq](https://github.com/pachterlab/ffq) to download the relevant data from ReMap. 
