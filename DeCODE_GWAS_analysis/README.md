# DeCODE GWAS Analysis pipeline
This pipeline carries out the standard preliminary analysis of GWAS results.
It begins from the files recieved from DeCODE Genetics as they are, re-formats them and extracts information from them.

The copying and unzipping of the files is taken care of by Snakemake. See more information about the pipeline in the [Snakefile](Snakefile), and more information about Snakemake below.
## Scripts
### **extract_variants_by_pval.py**
This script takes both the ´variant_info.txt´ file and and the chi2 tables and filters out entries based on a significance threshold.
### **produce_full_sumstats_for_single_trait.py**

## About Snakemake
[Snakemake](https://snakemake.github.io/) is a workflow management system that ensures the robustness and reproducibility of the analysis.
In order to run the pipeline in your computer, you need to have Snakemake installed. You can install it using conda.
You will also need to update the ´config.yaml´ file with the paths to the files you want to use.