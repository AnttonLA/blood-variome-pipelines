def check_index(bgens_list: list) -> None:
    """
    Check if the BGEN files are indexed. If they are not, ask the user if they want to index them.

    :param bgens_list: List of absolute paths of the BGEN files
    """
    lacking_index = False
    for bgen in bgens_list:
        if not os.path.isdir(bgen + '.idx2'):
            lacking_index = True
            break

    # If any of the bgen files are not indexed, ask the user if they want to index them
    if lacking_index:
        print("At least some of the bgen files are not indexed. Do you want to index them? (y/n)")
        while True:
            answer = input()
            if answer == 'y':
                hl.index_bgen(bgens_list)
                break
            elif answer == 'n':
                print("Exiting")
                sys.exit(0)
            else:
                print("Please answer 'y' or 'n'")
    else:
        print("All bgen files are indexed. Proceeding...")


def sanitize_filename(filename: str) -> str:
    """
    Remove invalid characters from a filename. Used when we want to use a column from a dataframe as a filename to make
    sure we will get a valid name.
    """
    invalid_chars = r' <>:"/\|?*'
    valid_chars = [c for c in filename if c not in invalid_chars]
    sanitized_filename = ''.join(valid_chars).strip(' .')
    return sanitized_filename


def extract_chromosome_number(file_path: str):
    """Extract the chromosome number from a file path."""
    # Use regular expression to find the number after 'chr' and before '.bgen'
    match = re.search(r'chr(\d+)\.bgen', file_path)
    if match:
        return int(match.group(1))
    else:
        return None


def sort_bgen_files_by_chromosome(bgen_paths: list) -> list:
    """Sort a list of BGEN file paths by chromosome number."""
    return sorted(bgen_paths, key=extract_chromosome_number)
