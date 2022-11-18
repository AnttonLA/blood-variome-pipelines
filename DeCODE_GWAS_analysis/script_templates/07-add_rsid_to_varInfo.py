import polars as pl
import tabix


dbsnp_path = "/home/antton/cbio3/data/dbSNP/GCF_000001405.39.gz"
path_to_chr_RefSeq = '/'.join(
    dbsnp_path.split('/')[:-1] + ['chr_to_RefSeq.txt'])  # Path to file with chr to RefSeq mappings
with open(path_to_chr_RefSeq, 'r') as f:  # Load the file with the mappings into a dictionary
    chr_to_RefSeq_dict = {line.split('\t')[0]: line.split('\t')[1].rstrip() for line in f}


def dbsnp_single_position_query(SNP_chr, SNP_pos: int):
    """
    Query the dbSNP file for a single position.
    It converts the chromosome number to the RefSeq chromosome name before querying.

    :param SNP_chr: chromosome number or name
    :param SNP_pos: position on the chromosome
    :return: full row from the dbSNP database corresponding to that position as a list
    """
    if SNP_chr == 'X':
        return None
    tb = tabix.open(dbsnp_path)
    query_str = f"{chr_to_RefSeq_dict[SNP_chr]}:{SNP_pos}-{SNP_pos}"  # For example: NC_000006.12:17100-17100
    matches = tb.querys(query_str)
    # The columns in the dbSNP file are: CHROM POS ID REF ALT QUAL FILTER INFO
    match_list = [x for x in matches]  # Convert the generator to a list
    if not match_list:
        print(f"No matches for {SNP_chr}:{SNP_pos}-{SNP_pos}")
        return None
    else:
        if len(match_list) > 1:
            for match in match_list:  # Loop through matches and take the one with position identical to SNP_pos
                match_position = int(match[1])
                if match_position == SNP_pos:  # Ensure the variant is in the exact same position that was queried
                    match_list = [match]  # Make it a "list of lists" even if it's only one list
                    break

        return match_list[0][2]  # return the rsID


df = pl.read_csv("variant_info_extended.txt", sep="\t")
#df = df.sample(40)
rsid_col = df.apply(lambda x: dbsnp_single_position_query(x[7][3:], x[8]))
df.insert_at_idx(1, rsid_col.to_series())
df = df.rename({'apply': 'rsid'})
print("DONE!")
print(df.head())
df.write_csv("variant_info_extended_with_rsid.txt", sep="\t", has_header=True)
