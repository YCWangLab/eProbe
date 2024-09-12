import pandas as pd


def split_lists(lists):
    first_elements = []
    second_elements = []
    for lst in lists:
        try:
            first_elements.append(lst[0])
            second_elements.append(lst[1])
        except IndexError:
            raise IndexError('''List index out of range, Use ', ' to separate upper and lower limits.''')

    return first_elements, second_elements


def filter_df_by_range(df, column_names, lower_limits, upper_limits):
    if len(column_names) != len(lower_limits) or len(column_names) != len(upper_limits):
        raise ValueError('Lengths of column_names, lower_limits, and upper_limits must be equal.')

    filtered_df = df.copy()
    for column_name, lower_limit, upper_limit in zip(column_names, lower_limits, upper_limits):
        if column_name in filtered_df.columns:
            print(f'''Filtering '{column_name}' based on range with limits ({lower_limit}, {upper_limit})''')
            filtered_df = filtered_df[
                (filtered_df[column_name] >= lower_limit) & (filtered_df[column_name] <= upper_limit)
                ]
    return filtered_df


def filter_df_by_percentile(df, column_name, lower_limit, upper_limit):
    if not 0 <= lower_limit <= upper_limit <= 1:
        raise ValueError('Percentage limit must be between 0 and 1. Upper limit must lager than lower limit.')

    print(f'''Filtering '{column_name}' based on percentile with limits ({lower_limit}, {upper_limit})''')

    sorted_df = df.sort_values(by=column_name)

    total_rows = sorted_df.shape[0]
    sorted_df['Percentile'] = (sorted_df.reset_index().index + 1) / total_rows

    filter_df = sorted_df[(lower_limit <= sorted_df['Percentile']) & (sorted_df['Percentile'] <= upper_limit)]
    filter_df = filter_df.drop(columns=['Percentile'])

    return filter_df


def parse_range_arg(arg):
    return list(map(float, arg.split(',')))


if __name__ == '__main__':
    import argparse

    usage = ''' python Tags_filter.py -m merged_tags.tsv --tm 60,80 --gc 40,60 
    --hairpin 0,30 --dust 0,2 --dimer 0,0.85 -o output_prefix'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-m', '--merged_tags', dest='merged_tags', action='store', nargs='?',
                        help='Merged tagged SNPs dataframe.',
                        metavar='FILE', required=True)
    parser.add_argument('-o', '--out', dest='output', action='store', nargs='?',
                        help='Filtered tagged SNPs data frame.',
                        metavar='STRING', required=True)
    parser.add_argument('--tm', dest='keep_tm', action='store', nargs='?',
                        help='Upper and lower limits for melting temperature (default: 68,78).',
                        default='60,80', metavar='STRING')
    parser.add_argument('--gc', dest='keep_gc', action='store', nargs='?',
                        help='Upper and lower limits for GC content (default: 40,60).',
                        default='40,60', metavar='STRING')
    parser.add_argument('--hairpin', dest='keep_hairpin', action='store', nargs='?',
                        help='Upper and lower limits for hairpin score (default: 0,60).',
                        default='0,60', metavar='STRING')
    parser.add_argument('--complexity', dest='keep_complexity', action='store', nargs='?',
                        help='Upper and lower limits for S score (sequence complexity) (default: 0,2).',
                        default='0,2', metavar='STRING')
    parser.add_argument('--dimer', dest='keep_dimer', action='store', nargs='?',
                        help='Upper and lower percentile limits for dimer structure (default: 0,0.85)',
                        default='0,0.85', metavar='STRING')
    args = parser.parse_args()

    snp_df = pd.read_csv(args.merged_tags, delimiter='\t', header=0)

    keep_tm, keep_gc, keep_hairpin, keep_dust = [parse_range_arg(args.keep_tm), parse_range_arg(args.keep_gc),
                                                 parse_range_arg(args.keep_hairpin),
                                                 parse_range_arg(args.keep_complexity)]
    label_names = ['tm', 'gc', 'hairpin', 'complexity']
    limits_list = [keep_tm, keep_gc, keep_hairpin, keep_dust]
    lower_limits, upper_limits = split_lists(limits_list)
    filtered_df = filter_df_by_range(snp_df, label_names, lower_limits, upper_limits)

    if 'dimer' in filtered_df.columns:
        keep_dimer = parse_range_arg(args.keep_dimer)
        filtered_df = filter_df_by_percentile(filtered_df, 'dimer', *keep_dimer)

    filtered_df.sort_index().to_csv(args.output + '.tsv', sep='\t', index=False)
    remove_SNP = int(snp_df.shape[0]) - (filtered_df.shape[0])
    print(f'Filtered out {remove_SNP} from {snp_df.shape[0]} probes in biophysical/biochemical filtering.')
