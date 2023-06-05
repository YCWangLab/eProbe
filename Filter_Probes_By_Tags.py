import pandas as pd
import argparse

usage = """python Filter_Probes_By_Tags.py -p pre_probes.fasta -g GCTm.txt -a hairpin.txt -d dimer.txt -u dust.txt -t 0 -o filtered_tags.txt"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-m", "--merged_labels", dest="merged_labels", action="store", nargs='?',
                    help="Merged labels file (.txt)", metavar="FILE")
parser.add_argument("-o", "--out", dest="output", action="store", nargs='?',
                    help="SNPs dataframe with filtered labels file (.txt)", metavar="FILE")
parser.add_argument("-t", "--Tm_limit", dest="keep_Tm", action="store", nargs='?',
                    help="Upper and lower limits for melting temperature (68,78)", metavar="FILE")
parser.add_argument("-g", "--GC_limit", dest="keep_GC", action="store", nargs='?',
                    help="Upper and lower limits for GC content (40,60)", metavar="FILE")
parser.add_argument("-a", "--hairpin_limit", dest="keep_hairpin", action="store", nargs='?',
                    help="Upper and lower limits for hairpin structure (0,30)", metavar="FILE")
parser.add_argument("-u", "--dust_limit", dest="keep_dust", action="store", nargs='?',
                    help="Upper and lower limits for complexity (0,2)", metavar="FILE")
parser.add_argument("-d", "--dimer_ratio", dest="keep_dimer", action="store", nargs='?',
                    help="Lower limits for dimer structure (0.85)", metavar="FILE")
parser.add_argument("-i", "--keep_taxid", dest="keep_taxid", action="store", nargs='?',
                    help="Taxid to be kept (0,4050,4051; 0 means failing to be classified by Kraken2)", metavar="FILE")
args = parser.parse_args()
# parameters
merged_labels = args.merged_labels
output_file = args.output
# keep Tm
keep_Tm = list(map(float, args.keep_Tm.split(",")))
# keep GC
keep_GC = list(map(float, args.keep_GC.split(",")))
# keep hairpin
keep_hairpin = list(map(float, args.keep_hairpin.split(",")))
# keep DUST score
keep_dust = list(map(float, args.keep_dust.split(",")))
# keep dimer upper limit ratio (0-1)
keep_dimer = float(args.keep_dimer)
# keep taxid
keep_taxid = args.keep_taxid.split(",")


# make range lists
def split_lists(lists):
    first_elements = []
    second_elements = []
    for lst in lists:
        first_elements.append(lst[0])
        second_elements.append(lst[1])

    return first_elements, second_elements


# keep probes by range
def keep_by_range(df, column_names, lower_limits, upper_limits):
    filtered_df = df.copy()
    for column_name, lower_limit, upper_limit in zip(column_names, lower_limits, upper_limits):
        filtered_df = filtered_df[
            (filtered_df[column_name] >= lower_limit) & (filtered_df[column_name] <= upper_limit)
            ]
    return filtered_df


# keep probes by percentage
def keep_by_percent(df, column_name, percentage_limit):
    sorted_df = df.sort_values(by=column_name)
    # ranking
    total_rows = sorted_df.shape[0]
    sorted_df['Percentile'] = (sorted_df.reset_index().index + 1) / total_rows
    # filter
    filtered_df = sorted_df[sorted_df['Percentile'] <= percentage_limit]
    filtered_df.drop(columns=['Percentile'], inplace=True)
    return filtered_df


# keep probes by values
def keep_by_values(df, column_name, target_list):
    # target column to list
    df_copy = df.copy()
    df_copy[column_name] = df_copy[column_name].str.split(',')
    # if subset or equal
    mask = df_copy[column_name].apply(lambda x: set(x) <= set(target_list))

    # return bool list
    filtered_df = df[mask]

    return filtered_df


# read dataframe
merged_df = pd.read_csv(merged_labels, delimiter='\t', header=0)

# filter Tm GC hairpin DUST by range
label_names = ["Tm", "GC", "hairpin", "DUST"]
limits_list = [keep_Tm, keep_GC, keep_hairpin, keep_dust]
lower_limits, upper_limits = split_lists(limits_list)
filtered_df1 = keep_by_range(merged_df, label_names, lower_limits, upper_limits)

# filter dimer by percentage
filtered_df2 = keep_by_percent(filtered_df1, 'dimer', keep_dimer)

# filter taxid with values
filtered_df3 = keep_by_values(filtered_df2, 'taxid', keep_taxid)

# output
df_sort = filtered_df3.sort_index()
df_sort.to_csv(output_file, sep='\t', index=False)
