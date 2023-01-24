# DeCODE GWAS Analysis pipeline
This pipeline carries out the standard preliminary analysis of GWAS results.
It begins from the files received from DeCODE Genetics as they are, re-formats them and extracts information from them.

The copying and unzipping of the files is taken care of by Snakemake. See more information about the pipeline in the
[Snakefile](Snakefile), and more information about Snakemake below.

## Usage
Each of the scripts in this pipeline can be run separately. Most of them use
[argparse](https://docs.python.org/3/library/argparse.html) to parse the arguments, so you can see the arguments by
running the script with the -h flag.

However, it is recommended that you run the full pipeline as one. You can accomplish this using
[Snakemake](https://snakemake.github.io/).

### About Snakemake
Snakemake is a workflow management system that ensures the robustness and reproducibility of the analysis.
In order to run the pipeline in your computer, you need to have Snakemake installed. You can install it through
[Conda](https://docs.conda.io/en/latest/) with the following commands:

    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install snakemake

You will also need to update the [config.yaml](config.yaml) file with the paths to the files you want to use.

Once you have updated the congig.yaml file so that it points to the correct files and directories, you can run the
pipeline with the following command:

    snakemake --cores 4

## Scripts
### **extract_variants_by_pval.py**
This script takes both the ´variant_info.txt´ file and the chi2 tables and filters out entries based on a significance
threshold.

### **produce_full_sumstats_for_single_trait.py**
This script takes the ´variant_info.txt´ file and the chi2 tables and produces a complete summary statistics file for a
single trait.

### **produce_all_hits_table.py**
This script takes the ´variant_info.txt´ file and the chi2 tables and produces a table with all hits for a single trait.

**Note regarding "repeats":** The default option is to NOT take repeats, meaning the script will take only the trait
that is most significant for each unique position/ID. In this version of the output table of hits, each row belongs to a
distinct variant, showing the stats of the phenotype that had the lowest p-value for that SNP. If the option to keep the
repeats is selected instead, there will be a separate entry for the same variant for each of the different phenotypes
where the variant has a significant association.  

### **swap_in_hits_into_template_sumstats_file.py**
This script takes a template summary statistics file and a table of hits and swaps in the hits into the template file.

### **plot_manhattan.sh**
This script takes a summary statistics file and plots a Manhattan plot for it using
[manhattan_maker](https://github.com/AnttonLA/manhattan_maker). You will need manhattan_maker installed on your
computer to use it.

### **generate_hit_regions_bed.py**
This script is used to generate a .bed file that contains the regions where at least one GWAS hit is present.
