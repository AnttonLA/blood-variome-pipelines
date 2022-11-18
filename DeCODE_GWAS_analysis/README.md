# DeCODE GWAS Analysis pipeline
This pipeline carries out the standard preliminary analysis of GWAS results.
It begins from the files recieved from DeCODE Genetics as they are, re-formats them and extracts information from them.

### **extract_variants_by_pval.py**
This script takes both the ´variant_info.txt´ file and and the chi2 tables and filters out entries based on a significance threshold.

