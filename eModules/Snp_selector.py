import collections
import pandas as pd
import multiprocessing
import random
from pipe_controler import check_col_name_in_df


def sliding_window_chr(chr_len, win_size):
    start = 1
    end = win_size
    result = []

    while end < chr_len:
        result.append([start, end])
        start += win_size
        end += win_size

    result.append([start, min(chr_len, end)])

    return result


def sliding_window_genome(genome_size_file, win_size):
    win_dict = collections.OrderedDict()
    with open(genome_size_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                chr_id, length = parts[0], int(parts[1])
                win_dict[chr_id] = sliding_window_chr(length, win_size)
    return win_dict


def split_dataframe(df, num_parts):
    if num_parts <= 0 or num_parts > len(df):
        raise ValueError('Invalid number of parts.')

    num_rows = len(df)
    rows_per_part = num_rows // num_parts

    df_parts = []
    start = 0

    for _ in range(num_parts - 1):
        end = start + rows_per_part
        df_parts.append(df.iloc[start:end])
        start = end

    df_parts.append(df.iloc[start:])
    return df_parts


def split_dict_equally(input_dict, num_parts):
    it = iter(input_dict)
    split_dicts = [{} for _ in range(num_parts)]
    for i, key in enumerate(it):
        split_dicts[i % num_parts][key] = input_dict[key]

    return split_dicts


def read_bed_df(bed_file_):
    return pd.read_table(bed_file_, sep='\t', header=None)


def keep_snp_by_sub_bed(sub_bed_df, snp_df, feature=None):
    sub_snp_keep_list = []
    for row in sub_bed_df.itertuples(index=False):
        chrom, start, end, region_type = row[0], int(row[1]), int(row[2]), row[3] if len(row) > 3 else None
        # specified feature and interval matched
        if feature and region_type == feature:
            snps_within_region = snp_df[(snp_df['chr'] == chrom) & (snp_df['pos'].between(start, end))]
            sub_snp_keep_list.extend(snps_within_region.index.tolist())
        # no specified feature, keep all within the interval
        elif not feature:
            snps_within_region = snp_df[(snp_df['chr'] == chrom) & (snp_df['pos'].between(start, end))]
            sub_snp_keep_list.extend(snps_within_region.index.tolist())
    return sub_snp_keep_list


def keep_snp_by_bed(snp_df_, bed_df_, thread, col_num, feature=None):
    if col_num >= 0:
        bed_df_.rename(columns={bed_df_.columns[col_num]: 'type'}, inplace=True)
    else:
        print('No column of BED was specified to filter by. Keeping all SNPs in intervals of BED.')

    if not feature:
        print('No feature was specified to filter by. Keeping all SNPs in intervals of BED.')

    all_snp_keep_list = []
    with multiprocessing.Pool(processes=int(thread)) as pool:
        split_dfs = split_dataframe(bed_df_, thread * 2)
        fixed_values = (snp_df_, feature)
        parameters_list = [(df,) + fixed_values for df in split_dfs]
        results = [pool.apply_async(keep_snp_by_sub_bed, args=params) for params in parameters_list]
        for result in results:
            all_snp_keep_list.extend(result.get())
        all_snp_keep_list = sorted(list(set(all_snp_keep_list)))
    return snp_df_.loc[all_snp_keep_list], all_snp_keep_list


def sort_df_by_columns(df, column_order):
    sorted_columns = [col for col in column_order if col in df.columns]
    return df.reindex(columns=sorted_columns)


def select_SNP_with_weight(df, weights, col_order):
    df_copy = df.copy()
    df_copy.drop(columns=['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt'], inplace=True)
    sorted_df = sort_df_by_columns(df_copy, col_order)
    normalized_df = (sorted_df - sorted_df.min()) / (sorted_df.max() - sorted_df.min())
    reversed_df = 1 - normalized_df
    weighted_df = reversed_df.multiply(weights)
    return weighted_df.sum(axis=1).idxmax()


def select_SNP_in_windows(window_dict, snp_df_, col_order, weights_=None, feature_index_list=None):
    selected_indexes = []
    for chr, windows in window_dict.items():
        for win_start, win_end in windows:
            win_SNPs_index_list = snp_df_[
                (snp_df_['chr'] == chr) & (snp_df_['pos'] >= win_start) & (snp_df_['pos'] <= win_end)
                ].index.tolist()

            # no SNPs in this window, skipping window
            if not win_SNPs_index_list:
                continue

            # has input feature BED
            if feature_index_list:
                intersec_list = list(set(win_SNPs_index_list) & set(feature_index_list))
                # has feature SNPs in this window, prioritizing feature SNPs
                if intersec_list:
                    selected_win_index = select_SNP_with_weight(snp_df_.loc[intersec_list], weights_,
                                                                col_order) if weights_ else random.choice(intersec_list)
                # no feature SNPs in this window, random selection or optimization in this window
                else:
                    selected_win_index = select_SNP_with_weight(snp_df_.loc[win_SNPs_index_list], weights_,
                                                                col_order) if weights_ else random.choice(
                        win_SNPs_index_list)
            # no input feature BED, random selection or optimization in all windows
            else:
                selected_win_index = select_SNP_with_weight(snp_df_.loc[win_SNPs_index_list], weights_,
                                                            col_order) if weights_ else random.choice(
                    win_SNPs_index_list)

            selected_indexes.append(selected_win_index)

    return snp_df_.loc[sorted(selected_indexes)]


def process_sub_windows(args):
    return select_SNP_in_windows(*args)


if __name__ == '__main__':
    import argparse

    usage = '''python Snp_selector.py -d snp_df.tsv -t desire_tm -g desire_gc -w weights -c chr_len.txt -s win_size 
    -b bed.txt -n column_num -f feature -m mode -r thread -o output_prefix'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-d', '--snp_df', dest='snp_df', action='store', nargs='?',
                        help='Filtered tags file.', metavar='FILE', required=True)
    parser.add_argument('-t', '--desire_tm', dest='desire_Tm', action='store', nargs='?',
                        help='Desire melting temperature (default: 75).',
                        default='75', metavar='STRING')
    parser.add_argument('-g', '--desire_gc', dest='desire_GC', action='store', nargs='?',
                        help='Desire GC content (default: 50).',
                        default='50', metavar='STRING')
    parser.add_argument('-w', '--weights', dest='weights', action='store', nargs='?',
                        help='Weights for gc, tm, complexity, hairpin, dimer (e.g., 0.5,0.5,0,0,0)',
                        metavar='STRING')
    parser.add_argument('-c', '--g_size', dest='g_size', action='store', nargs='?',
                        help='Tab-separated genome size file indicating the id and length of chromosomes designed on.',
                        metavar='FILE', required=True)
    parser.add_argument('-s', '--win_size', dest='win_size', action='store', nargs='?', default='1000',
                        help='Window size (default: 1000)', metavar='STRING')
    parser.add_argument('-b', '--bed', dest='bed', action='store', nargs='?',
                        help='Annotation file (bed format, converted from gff).', metavar='FILE')
    parser.add_argument('-n', '--column', dest='col_num', action='store', nargs='?',
                        help='Number of which column in BED file to be used to match feature.',
                        default='0', metavar='STRING')
    parser.add_argument('-f', '--feature', dest='feature', action='store', nargs='?',
                        help='Which feature in the column want to be kept (e.g., exon).',
                        metavar='STRING')
    parser.add_argument('-m', '--mode', dest='mode', action='store', nargs='?',
                        help='Select mode [both, feature, window] (default: window).',
                        default='window', metavar='STRING')
    parser.add_argument('-r', '--thread', dest='thread', action='store', nargs='?', default=1,
                        help='Number of thread (default: 1).', metavar='STRING')
    parser.add_argument('-o', '--output_prefix', dest='output', action='store', nargs='?',
                        help='Prefix for output', metavar='STRING')
    args = parser.parse_args()

    sliding_win_dict = sliding_window_genome(args.g_size, int(args.win_size))
    tag_order = ['gc', 'tm', 'complexity', 'hairpin', 'dimer']
    snp_df = pd.read_csv(args.snp_df, delimiter='\t', header=0)

    if args.weights:
        snp_df['tm'] = (snp_df['tm'] - float(args.desire_Tm)).abs()
        snp_df['gc'] = (snp_df['gc'] - float(args.desire_GC)).abs()
        weights = list(map(float, args.weights.split(',')))
        num_nonzero_weights = sum(1 for num in weights if num != 0.0)
        num_tags = sum(1 for col in snp_df.columns if col in tag_order)

        if abs(sum(weights) - 1.0) > 1e-9:
            raise ValueError('Sum of weights is not equal to 1.')
        elif num_nonzero_weights > num_tags:
            raise ValueError(
                f'''Number of weights ({num_nonzero_weights}) don't match the number of tags ({num_tags}) in {args.snp_df}.''')
        else:
            tag_weights = ', '.join(
                [f'{tag_weight[0]}:{tag_weight[1]}' for tag_weight in list(zip(tag_order, weights))])
            print(f'Enable tag-optimizer in a window with weights {tag_weights}.')
    else:
        weights = None
        print(f'Turn off tag-optimizer in a window since weights unset. Skipping...')

    if not check_col_name_in_df(tag_order, snp_df):
        weights = None
        print(f'Turn off tag-optimizer in a window since no tags are detected in {args.snp_df}.')

    if args.bed:
        anno_bed_df = read_bed_df(args.bed)

        if args.col_num == '0' or not args.feature:
            print(
                'Since you haven\'t specified which column or feature within the BED to filter by, all SNPs within "Bed" will be retained.')

        feature_SNPs_df, feature_SNPs_indexes = keep_snp_by_bed(snp_df, anno_bed_df, int(args.thread),
                                                                int(args.col_num) - 1, args.feature)
    else:
        feature_SNPs_df, feature_SNPs_indexes = None, None

    if args.mode == 'feature' and feature_SNPs_df is not None and not feature_SNPs_df.empty:
        feature_SNPs_df.drop(columns=tag_order, inplace=True)
        feature_SNPs_df.to_csv(f'{args.output}_feature_SNPs.tsv', sep='\t', index=False)
    elif args.mode == 'window':
        window_dict_chunks = split_dict_equally(sliding_win_dict, num_parts=int(args.thread))
        sliding_tasks = [(chunk, snp_df, tag_order, weights, feature_SNPs_indexes) for chunk in window_dict_chunks]
        with multiprocessing.Pool(processes=int(args.thread)) as pool:
            results = pool.map(process_sub_windows, sliding_tasks)
        window_SNPs_df = pd.concat(results)
        window_SNPs_df.drop(columns=tag_order, inplace=True)
        window_SNPs_df.to_csv(f'{args.output}_window_SNPs.tsv', sep='\t', index=False)
    elif args.mode == 'both' and feature_SNPs_df is not None and not feature_SNPs_df.empty:
        window_dict_chunks = split_dict_equally(sliding_win_dict, num_parts=int(args.thread))
        sliding_tasks = [(chunk, snp_df, tag_order, weights, feature_SNPs_indexes) for chunk in window_dict_chunks]
        with multiprocessing.Pool(processes=int(args.thread)) as pool:
            results = pool.map(process_sub_windows, sliding_tasks)
        window_SNPs_df = pd.concat(results)
        merged_SNPs_df = pd.merge(feature_SNPs_df, window_SNPs_df, how='outer')
        sorted_merged_SNPs_df = merged_SNPs_df.sort_index()

        sorted_merged_SNPs_df.drop(columns=tag_order, inplace=True)
        sorted_merged_SNPs_df.to_csv(f'{args.output}_merged_SNPs.tsv', sep='\t', index=False)
        feature_SNPs_df.drop(columns=tag_order, inplace=True)
        feature_SNPs_df.to_csv(f'{args.output}_feature_SNPs.tsv', sep='\t', index=False)
        window_SNPs_df.drop(columns=tag_order, inplace=True)
        window_SNPs_df.to_csv(f'{args.output}_window_SNPs.tsv', sep='\t', index=False)
    else:
        raise ValueError('Input proper select mode: both/feature/window.')
