#!/bin/bash

# This script calls manhattan_maker to create a Manhattan plot for the data in the input sumstats file.
# If manhattan_maker is not installed, it will attempt to install it.
#
# Usage: plot_manhattan.sh <path to manhattan_maker> <input_file> <output_file>
#
# Example command:
# ./plot_manhattan.sh "/home/antton/Projects/Immune_GWAS/data/processed/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/combined_manhattan.txt" "/home/antton/Projects/Immune_GWAS/data/processed/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/output_manhattan"

# Check that all inputs are provided:
if [ $# -ne 2 ]; then
    echo "Incorrect number of inputs. Usage: plot_manhattan.sh <input_file> <output_file>"
    exit 1
fi

# Check if manhattan_maker is installed
if ! [ -x "$(command -v manhattan_maker)" ]; then
    echo "manhattan_maker is not installed. Please install manhattan_maker before running this script." >&2
    exit 1
else
    echo "manhattan_maker is installed. Continuing..."
fi


# $1 E.g. "/home/antton/Projects/Immune_GWAS/data/processed/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/combined_manhattan.txt"\
# $2 E.g. -o "/home/antton/Projects/Immune_GWAS/data/processed/BloodVariome_Taravero_preliminary_GWAS_2022-10-29/Frequency_and_Ratio/output_manhattan"

#manhattan_maker --in-file $1 --col-name "ID" --col-chr "chromosome" --col-pos "position" --col-pvalue "pval"\
# --significant-threshold 11 --point-size 4 --significant-point-size 4 --abline 7.301029995664,11\
# --odd-chromosome-color "#ececec" --even-chromosome-color "#a9a9a9" --significant-color "#808080"\
# --chromosome-box-color "#ffffff" --axis-text-size 13 --chr-text-size 13 --label-text-size 28 --graph-height 6\
# --graph-width 15 -o $2

 manhattan_maker --in-file $1 --col-name "ID" --col-chr "chromosome" --col-pos "position" --col-pvalue "pval"\
 --significant-threshold 11 --point-size 4 --significant-point-size 4 --abline 7.301029995664,11\
 --odd-chromosome-color "#ececec" --even-chromosome-color "#a9a9a9" --significant-color "#808080"\
 --use-groups --col-groups "alias" --group-colors "#E95D1F,#3C8D2E,#9C00E9,#00E91B,#E900E1,#00E0E9,#E9C800"\
 --chromosome-box-color "#ffffff" --axis-text-size 13 --chr-text-size 13 --label-text-size 28 --graph-height 6\
 --graph-width 15 -o $2
