import pandas as pd


def split_lists(lists):
    first_elements = []
    second_elements = []
    for lst in lists:
        first_elements.append(lst[0])
        second_elements.append(lst[1])

    return first_elements, second_elements


def filter_df_by_range(df, column_names, lower_limits, upper_limits):
    if len(column_names) != len(lower_limits) or len(column_names) != len(upper_limits):
        raise ValueError("Lengths of column_names, lower_limits, and upper_limits must be equal.")

    filtered_df = df.copy()
    for column_name, lower_limit, upper_limit in zip(column_names, lower_limits, upper_limits):
        filtered_df = filtered_df[
            (filtered_df[column_name] >= lower_limit) & (filtered_df[column_name] <= upper_limit)
            ]
    return filtered_df


def filter_df_by_percentile(df, column_name, lower_limit, upper_limit):
    if not 0 <= lower_limit <= upper_limit <= 1:
        raise ValueError("Percentage limit must be between 0 and 1. Upper limit must lager than lower limit.")

    sorted_df = df.sort_values(by=column_name)

    total_rows = sorted_df.shape[0]
    sorted_df['Percentile'] = (sorted_df.reset_index().index + 1) / total_rows

    filter_df = sorted_df[(lower_limit <= sorted_df['Percentile']) & (sorted_df['Percentile'] <= upper_limit)]
    filter_df = filter_df.drop(columns=['Percentile'])

    return filter_df


def parse_range_arg(arg):
    return list(map(float, arg.split(",")))


if __name__ == '__main__':
    import argparse

    usage = """ python Tags_filter.py -m merged_tags.tsv --Tm 60,80 --GC 40,60 --hairpin 0,30 --dust 0,2 --dimer 0,0.85 -o output_prefix"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-m", "--merged_tags", dest="merged_tags", action="store", nargs='?',
                        help="Merged tagged SNPs dataframe.", metavar="FILE", required=True)
    parser.add_argument("-o", "--out", dest="output", action="store", nargs='?',
                        help="Filtered tagged SNPs dataframe.", metavar="STRING", required=True)
    parser.add_argument("--Tm", dest="keep_Tm", action="store", nargs='?', default="60,80",
                        help="Upper and lower limits for melting temperature (default: 68,78).", metavar="STRING")
    parser.add_argument("--GC", dest="keep_GC", action="store", nargs='?', default="40,60",
                        help="Upper and lower limits for GC content (default: 40,60). ", metavar="STRING")
    parser.add_argument("--hairpin", dest="keep_hairpin", action="store", nargs='?', default="0,30",
                        help="Upper and lower limits for hairpin score (default: 0,30).", metavar="STRING")
    parser.add_argument("--dust", dest="keep_dust", action="store", nargs='?', default="0,2",
                        help="Upper and lower limits for DUST score (sequence complexity) (default: 0,2).",
                        metavar="STRING")
    parser.add_argument("--dimer", dest="keep_dimer", action="store", nargs='?', default="0,0.85",
                        help="Upper and lower percentile limits for dimer structure (default: 0,0.85)", metavar="STRING")
    args = parser.parse_args()

    merged_df = pd.read_csv(args.merged_tags, delimiter='\t', header=0)

    keep_Tm, keep_GC, keep_hairpin, keep_dust = (parse_range_arg(args.keep_Tm), parse_range_arg(args.keep_GC),
                                                 parse_range_arg(args.keep_hairpin), parse_range_arg(args.keep_dust))
    label_names = ["Tm", "GC", "hairpin", "DUST"]
    limits_list = [keep_Tm, keep_GC, keep_hairpin, keep_dust]
    lower_limits, upper_limits = split_lists(limits_list)
    filtered_df = filter_df_by_range(merged_df, label_names, lower_limits, upper_limits)

    keep_dimer = parse_range_arg(args.keep_dimer)
    filtered_df = filter_df_by_percentile(filtered_df, 'dimer', *keep_dimer)

    filtered_df.sort_index().to_csv(args.output + ".tsv", sep='\t', index=False)
