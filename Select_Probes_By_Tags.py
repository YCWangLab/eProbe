import collections
import pandas as pd
import argparse
import multiprocessing

usage = """python Select_Probes_By_Tags.py -l filtered_tags.txt -t desire_Tm -g desire_GC -w weights -c chr_len.txt -s win_size -b bed.txt -n column_num -f feature -m mode -r threads -o output_prefix"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-l", "--filtered_tags", dest="filtered_tags", action="store", nargs='?',
                    help="Filtered labels file", metavar="FILE")
parser.add_argument("-t", "--desire_Tm", dest="desire_Tm", action="store", nargs='?',
                    help="Desire melting temperature", metavar="FILE")
parser.add_argument("-g", "--desire_GC", dest="desire_GC", action="store", nargs='?',
                    help="Desire GC content", metavar="FILE")
parser.add_argument("-w", "--weights", dest="weights", action="store", nargs='?',
                    help="Weights for Tm, GC, hairpin, dimer, dust (e.g., 0.8,0.2,0,0,0)", metavar="FILE")
parser.add_argument("-c", "--chr_length", dest="chr_len", action="store", nargs='?',
                    help="File storing length of each chromosome", metavar="FILE")
parser.add_argument("-s", "--win_size", dest="win_size", action="store", nargs='?',
                    help="Window size", metavar="FILE")
parser.add_argument("-b", "--bed", dest="bed", action="store", nargs='?',
                    help="Annotation file (bed format, converted from gff)", metavar="FILE")
parser.add_argument("-n", "--column", dest="col_num", action="store", nargs='?',
                    help="Which column in annotation file to be used to match feature", metavar="FILE")
parser.add_argument("-f", "--feature", dest="feature", action="store", nargs='?',
                    help="Which feature want to be kept (e.g., exon)", metavar="FILE")
parser.add_argument("-m", "--mode", dest="mode", action="store", nargs='?',
                    help="Selecting mode (both/feature/window)", metavar="FILE")
parser.add_argument("-r", "--threads", dest="threads", action="store", nargs='?',
                    help="Threads for both/feature mode", metavar="FILE")
parser.add_argument("-o", "--output_prefix", dest="output", action="store", nargs='?',
                    help="Prefix for output", metavar="FILE")
args = parser.parse_args()

# labels to be considered in this step: Tm(main) GC hairpin dimer dust
# parameters
filtered_labels = args.filtered_tags
output = args.output
# desire Tm
desire_Tm = float(args.desire_Tm)
# desire GC
desire_GC = float(args.desire_GC)
# weight for each labels (Tm, GC, hairpin, dimer, dust)
weights = list(map(float, args.weights.split(",")))
# chr length
chr_file = args.chr_len
# window size
window_size = int(args.win_size)
# annotation
anno_bed = args.bed
# col_num
col_num = int(args.col_num)
# feature
feature = args.feature
# mode (both/feature/window)
mode = args.mode
# threads
threads = int(args.threads)


def sliding_window_with_length(length, win_size):
    start = 1
    end = win_size
    result = []

    while end < length:
        result.append([start, end])
        start += win_size
        end += win_size
    result.append([start, length])
    return result


# get sliding windows list for whole genome
def sliding_window_with_txt(chrs_length, win_size=None):
    win_dict = collections.OrderedDict()
    with open(chrs_length, "r") as f:
        for line in f.readlines():
            chr = line.split()[0]
            length = line.split()[1]
            win_dict[chr] = sliding_window_with_length(int(length), win_size)
    f.close()
    return win_dict


# function for splitting dataframe
def split_dataframe(df, num_parts):
    num_rows = len(df)
    rows_per_part = num_rows // num_parts

    df_parts = []
    start = 0
    for i in range(num_parts - 1):
        end = start + rows_per_part
        df_part = df.iloc[start:end]
        df_parts.append(df_part)
        start = end

    df_part = df.iloc[start:]
    df_parts.append(df_part)

    return df_parts


# function for making annotation bed (converted from gff) df
def make_bed_df(bed_file_):
    # one-based coordinate
    bed_df = pd.read_table(bed_file_, sep="\t", header=None)
    return bed_df


def keep_snp_by_sub_bed(sub_bed_df, snp_df_, feature):
    sub_snp_keep_list = []
    for row in sub_bed_df.itertuples():
        chrom = row._1
        start = int(row._2)
        end = int(row._3)
        type_ = row.type
        # index list of snp within blocks
        if type_ == feature:
            bed_keep_list = snp_df_[
                (snp_df_.chr == chrom) & (snp_df_.pos >= start) & (snp_df_.pos <= end)].index.tolist()
            if len(bed_keep_list):
                sub_snp_keep_list = sub_snp_keep_list + bed_keep_list
    return sub_snp_keep_list


# function for keeping snp within given genome feature
def keep_snp_by_bed(snp_df_, bed_df_, col_num, feature, Thread):
    bed_df_.rename(columns={col_num: 'type'}, inplace=True)
    if __name__ == '__main__':
        with multiprocessing.Pool(processes=Thread) as pool:
            splitted_dfs = split_dataframe(bed_df_, Thread * 2)
            fixed_values = [snp_df_, feature]
            parameters_list = [[x] + fixed_values for x in splitted_dfs]
            All_snp_keep_list = []
            results = [pool.apply_async(keep_snp_by_sub_bed, args=(parameters)) for parameters in parameters_list]
            for result in results:
                All_snp_keep_list = All_snp_keep_list + result.get()
            All_snp_keep_list = sorted(list(set(All_snp_keep_list)))
            keep_snp_df = snp_df_.loc[All_snp_keep_list]
    return keep_snp_df, All_snp_keep_list


# select probes with weight labels
def select_row_with_weight(df, weights):
    df_copy = df.copy()
    df_copy.drop(columns=['chr', 'pos', 'type', 'taxid'], inplace=True)
    # Min-Max normalization
    normalized_df = (df_copy - df_copy.min()) / (df_copy.max() - df_copy.min())
    reversed_df = abs(normalized_df.sub(1))
    # weight
    weighted_df = reversed_df.multiply(weights)
    # select row
    selected_row = weighted_df.sum(axis=1).idxmax()
    return selected_row


# select probes in windows
def select_SNP_in_windows(window_dict, snp_df_, weights_, feature_index_list=None):
    selected_indexes = []
    for chr in window_dict.keys():
        for win in window_dict[chr]:
            win_start = int(win[0])
            win_end = int(win[1])
            # SNPs in windows
            win_SNPs_index_list = snp_df_[
                (snp_df_.chr == chr) & (snp_df_.pos >= win_start) & (snp_df_.pos <= win_end)].index.tolist()
            # intersection with exon SNPs, selecting from exon SNPs
            if feature_index_list and len(set(win_SNPs_index_list).intersection(feature_index_list)):
                intersec_list = list(set(win_SNPs_index_list).intersection(feature_index_list))
                selected_win_index = select_row_with_weight(snp_df_.loc[intersec_list], weights_)
            # without intersection or mode = window:
            else:
                if len(win_SNPs_index_list):
                    win_df = snp_df_.loc[win_SNPs_index_list]
                    selected_win_index = select_row_with_weight(win_df, weights_)
                else:
                    continue
            selected_indexes.append(selected_win_index)
    return snp_df_.loc[sorted(selected_indexes)]


# 1 reading filtered labels
filtered_labels_df = pd.read_csv(filtered_labels, delimiter='\t', header=0)
# distance to desire Tm
filtered_labels_df['Tm'] = (filtered_labels_df['Tm'] - desire_Tm).abs()
# distance to desire GC
filtered_labels_df['GC'] = (filtered_labels_df['GC'] - desire_GC).abs()
if mode == "both":
    # 2 making windows dict
    sliding_win_dict = sliding_window_with_txt(chr_file, window_size)
    # 3 reading bed file
    anno_bed_df = make_bed_df(anno_bed)
    # 4 All SNPs in Exon
    feature_SNPs_df, feature_SNPs_indexes = keep_snp_by_bed(filtered_labels_df, anno_bed_df, col_num - 1, feature, threads)
    # 5 SNPs in windows
    window_SNPs_df = select_SNP_in_windows(sliding_win_dict, filtered_labels_df, weights, feature_SNPs_indexes)
    # 6 merge
    merged_SNPs_df = pd.merge(feature_SNPs_df, window_SNPs_df, how='outer')
    sorted_merged_SNPs_df = merged_SNPs_df.sort_index()
    # output
    sorted_merged_SNPs_df.to_csv(output + "_merged_SNPs.txt", sep='\t', index=False)
    feature_SNPs_df.to_csv(output + "_%s_SNPs.txt" % feature, sep='\t', index=False)
    window_SNPs_df.to_csv(output + "_window_SNPs.txt", sep='\t', index=False)
elif mode == "feature":
    # 2 reading bed file
    anno_bed_df = make_bed_df(anno_bed)
    # 3 All SNPs in Exon
    exon_SNPs_df, exon_SNPs_indexes = keep_snp_by_bed(filtered_labels_df, anno_bed_df, col_num - 1, feature)
    exon_SNPs_df.to_csv(output + "_%s_SNPs.txt" % feature, sep='\t', index=False)
elif mode == "window":
    # 2 making windows dict
    sliding_win_dict = sliding_window_with_txt(chr_file, window_size)
    # 3 SNPs in windows
    window_SNPs_df = select_SNP_in_windows(sliding_win_dict, filtered_labels_df, weights)
    window_SNPs_df.to_csv(output + "_window_SNPs.txt", sep='\t', index=False)
