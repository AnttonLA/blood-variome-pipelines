#!/bin/bash

# This script call manhattan_maker.py to create Manhattan plots for all the
# samples in the combined_manhattan.txt file.

python manhattan_maker.py --in-file "/home/antton/Projects/Immune_GWAS/output/GWAS/combined_output/output_manhattan.txt"\
 --col-name "rsid" --col-chr "chromosome" --col-pos "position" --col-pvalue "pval" --significant-threshold 11\
 --point-size 4 --significant-point-size 4 --abline 7.301029995664,11 --odd-chromosome-color "#ececec"\
 --even-chromosome-color "#a9a9a9" --significant-color "#808080" --use-groups --col-group "lineage"\
 --group-colors "#b4a7d6,#80974d,#cd4646,#16367e,#232221" --chromosome-box-color "#ffffff" --axis-text-size 13\
 --chr-text-size 13 --label-text-size 28 --graph-height 6 --graph-width 15\
 -o "/home/antton/Projects/Immune_GWAS/output/GWAS/combined_output/output_manhattan_groups_test_4"
