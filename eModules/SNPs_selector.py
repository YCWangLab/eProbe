import collections
import pandas as pd
import multiprocessing
import random


# function: Generate sliding windows for a chromosome given its length and the window size
def sliding_window_chr(chr_len, win_size):
    """
    Generate sliding windows for a chromosome given its length and the window size.

    Parameters:
    - chr_len (int): Length of the chromosome.
    - win_size (int): Size of the sliding window.

    Returns:
    - list of lists: Each sub-list contains the start and end positions of a window.
    """
    start = 1
    end = win_size
    result = []

    # Generate windows until the end of the chromosome is reached
    while end < chr_len:
        result.append([start, end])
        start += win_size
        end += win_size

    # Adjusts the last window if smaller than the win_size
    result.append([start, min(chr_len, end)])

    return result


# function: Generates sliding windows for each chromosome in a genome
def sliding_window_genome(genome_size_file, win_size):
    """
    Generates sliding windows for each chromosome in a genome.

    Parameters:
    genome_size_file (str): The file path containing chromosome names and their lengths.
    win_size (int): The size of each sliding window.

    Returns:
    collections.OrderedDict: A dictionary with chromosome names as keys and lists of window positions as values.
    """
    win_dict = collections.OrderedDict()
    with open(genome_size_file, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                chr_id, length = parts[0], int(parts[1])
                # Generate windows for each chromosome
                win_dict[chr_id] = sliding_window_chr(length, win_size)
    return win_dict


# function: Splits a DataFrame into a specified number of parts as evenly as possible
def split_dataframe(df, num_parts):
    """
    Splits a DataFrame into a specified number of parts as evenly as possible.

    Parameters:
    - df (pandas.DataFrame): The DataFrame to split.
    - num_parts (int): The number of parts to divide the DataFrame into.

    Returns:
    - list of pandas.DataFrame: A list containing the split parts of the DataFrame.
    """
    if num_parts <= 0 or num_parts > len(df):
        raise ValueError("Invalid number of parts.")

    num_rows = len(df)
    rows_per_part = num_rows // num_parts

    df_parts = []
    start = 0

    # Split DataFrame into equal parts, except possibly the last one
    for _ in range(num_parts - 1):
        end = start + rows_per_part
        df_parts.append(df.iloc[start:end])
        start = end

    # Include any remaining rows in the last part
    df_parts.append(df.iloc[start:])
    return df_parts


# function: Split a dictionary into specified number of parts as evenly as possible
def split_dict_equally(input_dict, num_parts):
    """
    Split a dictionary into specified number of parts as evenly as possible.

    Parameters:
    - input_dict (dict): Dictionary to be split.
    - num_parts (int): Number of parts to split the dictionary into.

    Returns:
    - list of dicts: List containing the split parts of the dictionary.
    """
    it = iter(input_dict)
    split_dicts = [{} for _ in range(num_parts)]
    for i, key in enumerate(it):
        split_dicts[i % num_parts][key] = input_dict[key]

    return split_dicts


# function: Read annotation bed df (converted from gff)
def read_bed_df(bed_file_):
    # one-based coordinate
    bed_df = pd.read_table(bed_file_, sep="\t", header=None)
    return bed_df


# function: Filter SNPs within given genomic regions defined in a BED DataFrame and matching a specific feature
def keep_snp_by_sub_bed(sub_bed_df, snp_df, feature=None):
    """
    Filter SNPs within given genomic regions defined in a BED DataFrame and matching a specific feature.

    Parameters:
    - sub_bed_df (pd.DataFrame): DataFrame representing sub-parts of a BED file.
    - snp_df (pd.DataFrame): DataFrame containing SNP data.
    - feature (str): Feature type to match.

    Returns:
    - list: List of SNP indices that are within the given regions and match the feature.
    """
    sub_snp_keep_list = []
    for row in sub_bed_df.itertuples():
        chrom, start, end = row._1, int(row._2), int(row._3)
        snps_within_region = snp_df[(snp_df['chr'] == chrom) & (snp_df['pos'].between(start, end))]
        if feature and hasattr(snps_within_region, 'type'):
            snps_matching_feature = snps_within_region[snps_within_region['type'] == feature]
            sub_snp_keep_list.extend(snps_matching_feature.index.tolist())
        else:
            sub_snp_keep_list.extend(snps_within_region.index.tolist())
    return sub_snp_keep_list


# function: Filter SNPs within given genomic regions and matching a specific feature
def keep_snp_by_bed(snp_df_, bed_df_, thread, col_num, feature=None):
    """
    Filter SNPs within given genomic regions and matching a specific feature using multiprocessing.

    Parameters:
    - snp_df_ (pd.DataFrame): DataFrame containing SNP data.
    - bed_df_ (pd.DataFrame): DataFrame representing BED file.
    - col_num (int): Column number to match the feature.
    - feature (str): Feature type to match.
    - thread (int): Number of threads for multiprocessing.

    Returns:
    - pd.DataFrame: DataFrame of SNPs that are within the given regions and match the feature.
    """
    if col_num >= 0:
        bed_df_.rename(columns={col_num: 'type'}, inplace=True)
    all_snp_keep_list = []

    with multiprocessing.Pool(processes=int(thread)) as pool:
        # split the dataframe for multiprocessing
        split_dfs = split_dataframe(bed_df_, thread * 2)
        # make a parameter list for each sub_bed df
        fixed_values = (snp_df_, feature)
        parameters_list = [(df,) + fixed_values for df in split_dfs]
        results = [pool.apply_async(keep_snp_by_sub_bed, args=params) for params in parameters_list]
        for result in results:
            all_snp_keep_list.extend(result.get())
        all_snp_keep_list = sorted(list(set(all_snp_keep_list)))
        keep_snp_df = snp_df_.loc[all_snp_keep_list]
    return keep_snp_df, all_snp_keep_list


# function: Selects a single SNP from a DataFrame based on a weighted evaluation of its features
def select_SNP_with_weight(df, weights):
    """
    Selects a single SNP from a DataFrame based on a weighted evaluation of its features.
    Uses min-max normalization and weights to calculate a score for each SNP.

    Parameters:
    - df (pd.DataFrame): DataFrame of SNPs and their features.
    - weights (list): Weights to apply to each feature for scoring.

    Returns:
    - int: Index of the SNP with the highest score.
    """
    df_copy = df.copy()
    df_copy.drop(columns=["chr", "start", "end", "pos", "type", "ref", "alt"], inplace=True)
    # Min-Max normalization
    normalized_df = (df_copy - df_copy.min()) / (df_copy.max() - df_copy.min())
    # Invert normalized values
    reversed_df = 1 - normalized_df
    # Apply weights
    weighted_df = reversed_df.multiply(weights)
    # Select row with the highest sum of weighted values
    selected_row = weighted_df.sum(axis=1).idxmax()
    return selected_row


# function: Select a SNP from each window based on weighted features. Optionally filters SNPs to those intersecting with given features.
def select_SNP_in_windows(window_dict, snp_df_, weights_=None, feature_index_list=None):
    """
    Select a SNP from each window based on weighted features. Optionally filters SNPs to those intersecting with given features.

    Parameters:
    - window_dict (dict): Dictionary of windows per chromosome.
    - snp_df_ (pd.DataFrame): DataFrame of SNPs with 'chr' and 'pos' columns.
    - weights_ (list or pd.Series, optional): Weights for scoring SNPs.
    - feature_index_list (list, optional): List of SNP indices corresponding to a specific feature.

    Returns:
    - pd.DataFrame: DataFrame of selected SNPs.
    """
    selected_indexes = []
    for chr, windows in window_dict.items():
        for win_start, win_end in windows:
            # Identify SNPs within the current window
            win_SNPs_index_list = snp_df_[
                (snp_df_['chr'] == chr) & (snp_df_['pos'] >= win_start) & (snp_df_['pos'] <= win_end)
                ].index.tolist()

            # If no SNPs in the window, skip this iteration
            if not win_SNPs_index_list:
                continue

            # Determine the SNPs to consider for selection
            if feature_index_list:
                intersec_list = list(set(win_SNPs_index_list) & set(feature_index_list))
                if intersec_list:
                    if weights_:
                        selected_win_index = select_SNP_with_weight(snp_df_.loc[intersec_list],
                                                                                         weights_)
                    else:
                        # feature and random
                        selected_win_index = random.choice(intersec_list)
                else:
                    if weights_:
                        # weight
                        selected_win_index = select_SNP_with_weight(
                            snp_df_.loc[win_SNPs_index_list], weights_)
                    else:
                        # random
                        selected_win_index = random.choice(win_SNPs_index_list)
            else:
                if weights_:
                    # weight
                    selected_win_index = select_SNP_with_weight(snp_df_.loc[win_SNPs_index_list],
                                                                                     weights_)
                else:
                    # random
                    selected_win_index = random.choice(win_SNPs_index_list)
            selected_indexes.append(selected_win_index)

    # Return a DataFrame of the selected SNPs, sorted by index
    return snp_df_.loc[sorted(selected_indexes)]


# function: process a divided sub window list
def process_sub_windows(args):
    return select_SNP_in_windows(*args)


if __name__ == '__main__':
    import argparse

    usage = """python SNPs_selector.py -d snp_df.tsv -t desire_Tm -g desire_GC -w weights -c chr_len.txt -s win_size -b bed.txt -n column_num -f feature -m mode -r thread -o output_prefix"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-d", "--snp_df", dest="snp_df", action="store", nargs='?',
                        help="Filtered tags file.", metavar="FILE", required=True)
    parser.add_argument("-t", "--desire_Tm", dest="desire_Tm", action="store", nargs='?', default="75",
                        help="Desire melting temperature (default: 75).", metavar="STRING")
    parser.add_argument("-g", "--desire_GC", dest="desire_GC", action="store", nargs='?', default="50",
                        help="Desire GC content (default: 50).", metavar="STRING")
    parser.add_argument("-w", "--weights", dest="weights", action="store", nargs='?',
                        help="Weights for Tm, GC, hairpin, dimer, DUST (e.g., 0.5,0.5,0,0,0)", metavar="STRING")
    parser.add_argument("-c", "--g_size", dest="g_size", action="store", nargs='?',
                        help="Tab-separated genome size file indicating the id and length of chromosomes designed on.",
                        metavar="FILE", required=True)
    parser.add_argument("-s", "--win_size", dest="win_size", action="store", nargs='?', default="1000",
                        help="Window size (default: 1000)", metavar="STRING")
    parser.add_argument("-b", "--bed", dest="bed", action="store", nargs='?',
                        help="Annotation file (bed format, converted from gff).", metavar="FILE")
    parser.add_argument("-n", "--column", dest="col_num", action="store", nargs='?',
                        help="Number of which column in BED file to be used to match feature.",
                        default="0", metavar="STRING")
    parser.add_argument("-f", "--feature", dest="feature", action="store", nargs='?',
                        help="Which feature in the column want to be kept (e.g., exon).",
                        metavar="STRING")
    parser.add_argument("-m", "--mode", dest="mode", action="store", nargs='?',
                        help="Select mode [both, feature, window] (default: window).",
                        default="window",metavar="STRING")
    parser.add_argument("-r", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of thread (default: 1).", metavar="STRING")
    parser.add_argument("-o", "--output_prefix", dest="output", action="store", nargs='?',
                        help="Prefix for output", metavar="STRING")
    args = parser.parse_args()

    # make sliding windows dict
    sliding_win_dict = sliding_window_genome(args.g_size, int(args.win_size))

    # reading snp_df
    snp_df = pd.read_csv(args.snp_df, delimiter='\t', header=0)

    # weight for each tags (Tm, GC, hairpin, dimer, dust)
    if args.weights:
        # distance to desire Tm
        snp_df['Tm'] = (snp_df['Tm'] - float(args.desire_Tm)).abs()
        # distance to desire GC
        snp_df['GC'] = (snp_df['GC'] - float(args.desire_GC)).abs()
        weights = list(map(float, args.weights.split(",")))
    else:
        weights = None

    # bed for prioritarily keep SNP
    if args.bed:
        # reading annotation bed file
        anno_bed_df = read_bed_df(args.bed)
        # extract all SNPs within bed
        feature_SNPs_df, feature_SNPs_indexes = keep_snp_by_bed(snp_df, anno_bed_df, int(args.thread),
                                                                int(args.col_num) - 1, args.feature)
    else:
        feature_SNPs_df, feature_SNPs_indexes = None, None

    if args.mode == "feature" and feature_SNPs_df:
        feature_SNPs_df.to_csv(f"{args.output}_feature_SNPs.tsv", sep='\t', index=False)
    elif args.mode == "window":
        window_dict_chunks = split_dict_equally(sliding_win_dict, num_parts=int(args.thread))
        sliding_tasks = [(chunk, snp_df, weights, feature_SNPs_indexes) for chunk in window_dict_chunks]
        with multiprocessing.Pool(processes=int(args.thread)) as pool:
            results = pool.map(process_sub_windows, sliding_tasks)
        window_SNPs_df = pd.concat(results)
        window_SNPs_df.to_csv(f"{args.output}_window_SNPs.tsv", sep='\t', index=False)
    elif args.mode == "both" and feature_SNPs_df:
        window_dict_chunks = split_dict_equally(sliding_win_dict, num_parts=int(args.thread))
        sliding_tasks = [(chunk, snp_df, weights, feature_SNPs_indexes) for chunk in window_dict_chunks]
        with multiprocessing.Pool(processes=int(args.thread)) as pool:
            results = pool.map(process_sub_windows, sliding_tasks)
        window_SNPs_df = pd.concat(results)
        merged_SNPs_df = pd.merge(feature_SNPs_df, window_SNPs_df, how='outer')
        sorted_merged_SNPs_df = merged_SNPs_df.sort_index()
        sorted_merged_SNPs_df.to_csv(f"{args.output}_merged_SNPs.tsv", sep='\t', index=False)
        feature_SNPs_df.to_csv(f"{args.output}_feature_SNPs.tsv", sep='\t', index=False)
        window_SNPs_df.to_csv(f"{args.output}_window_SNPs.tsv", sep='\t', index=False)
