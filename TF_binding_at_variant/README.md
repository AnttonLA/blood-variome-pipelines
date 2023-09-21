# Transcription Factor Assignment Pipeline

Pipeline to find the Transcription Factor that is likely affected by a given non-coding genetic variant.

It takes as input a file with the following format:

```
ID	Chrom	Pos	OA	EA
variant_id	chr1	23935190	G	T
```

Currently, there are only two steps:

1. Look up entries of the ReMAP database (https://remap2022.univ-amu.fr/) to find Transcription Factors that
bind over a given variant.
2. Use PERFECTOS-APE (https://opera.autosome.org/perfectosape) to find the transcription factors whose binding motif likely disrupted by the variant.

The pipeline will then assign those TFs that appear in **BOTH** to the variant.

It **requires tabix** to perform the lookup on the ReMap metadata table and find all the relevant studies. 

### Usage

In order to run the pipeline, you will need:
* Input file with the format described above
* ReMap metadata BED file, bgzipped and tabix-indexed
* hg38 reference genome fasta file (used by PERFECTOS-APE)

You will also need to have **samtools installed**.

Edit the paths in the config.yaml file so that they point to the correct locations. 
Then, run the pipeline with:

```
snakemake --cores=all
```

### Future Additions

It would be nice to incorporate [ffq](https://github.com/pachterlab/ffq) to download the relevant data from ReMap. 

## References

E. Vorontsov, I.; V. Kulakovskiy, I.; Khimulya, G.; D. Nikolaeva, D. and J. Makeev, V. (2015). **PERFECTOS-APE - Predicting Regulatory Functional Effect of SNPs by Approximate P-value Estimation.** In *Proceedings of the International Conference on Bioinformatics Models, Methods and Algorithms (BIOSTEC 2015) - BIOINFORMATICS*; ISBN 978-989-758-070-3; ISSN 2184-4305, SciTePress, pages 102-108. DOI: 10.5220/0005189301020108
