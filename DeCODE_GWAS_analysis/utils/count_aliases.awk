# This script counts the number of different entries in the last column of a tab-separated table.
# It is meant to be used to calculate the number of unique aliases in the 'alias' column of the combined summary
# statistics file generated by the pipeline. Once we know this number, we can use it to dynamically generate a number
# of colors for plotting.
#
# Usage:
# awk -f count_aliases.awk combined_summary_stats.txt

# Set the input and output field separator to a tab character
BEGIN {
  FS = OFS = "\t"
}

# Store the entries in the last column in an array
{
  entries[$NF]++
}

# At the end, print the length of the array
END {
  print length(entries)
}