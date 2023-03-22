import polars as pl
import argparse


def convert_to_gor(input_file, chr_col, pos_col, additional_cols, sep, output_file):
    in_df = pl.read_csv(input_file, sep=sep, has_header=True)

    # Make sure the chr and pos column names exist in the input file
    if chr_col not in in_df.columns:
        raise Exception("Chromosome column not found: " + chr_col)
    if pos_col not in in_df.columns:
        raise Exception("Position column not found: " + pos_col)
    print("Input file:")
    print(in_df)
    # Reorder columns so that the order is chr, pos, additional columns
    cols = [chr_col, pos_col]
    if additional_cols:
        # Make sure the additional columns exist in the input file
        extra_col_list = []
        for col in additional_cols.split(","):
            col = col.strip()
            if col not in in_df.columns:
                raise Exception("Additional column not found: " + col)
            extra_col_list.append(col)
        cols.extend(extra_col_list)

    out_df = in_df.select(cols)
    # rename the columns to the GOR format: "#Chrom", "Pos"
    out_df = out_df.rename({chr_col: "#Chrom", pos_col: "Pos"})

    # Sort by chr and pos
    out_df = out_df.sort(["#Chrom", "Pos"])

    # check if '#Chrom' is a string column and convert it if not
    if not out_df['#Chrom'].dtype == pl.datatypes.Utf8:
        out_df = out_df.with_column(pl.col('#Chrom').cast(pl.datatypes.Utf8).alias("#Chrom"))
    # check if the elements of '#Chrom' start with 'chr'
    out_df = out_df.with_column(pl.when(~pl.col('#Chrom').str.starts_with('chr'))
            .then('chr' + pl.col('#Chrom'))
            .otherwise(pl.col('#Chrom')).alias("#Chrom"))

    print("Output file:")
    print(out_df)

    out_df.write_csv(output_file, sep="\t", has_header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input file")
    parser.add_argument("--chr_col", required=True, help="chromosome column", default="chr")
    parser.add_argument("--pos_col", required=True, help="position column", default="pos")
    parser.add_argument("--additional_cols", help="list of additional columns to include", default="")
    parser.add_argument("--sep", help="separator", default="\t")
    parser.add_argument("-o", "--output", help="output file")
    args = parser.parse_args()

    convert_to_gor(args.input, args.chr_col, args.pos_col, args.additional_cols, args.sep, args.output)
