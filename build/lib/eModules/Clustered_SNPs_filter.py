# -*- coding:utf-8 -*-
import pandas as pd
import numpy as np
from multiprocessing import Pool
import subprocess
import datetime
from pip_controler import diff_file_line_count

def count_surrounding_SNPs_chromosome(args):
    chromosome, positions, flank = args
    counts = []
    for pos in positions:
        distances = np.abs(positions - pos)
        count = np.sum(distances <= flank)
        counts.append(count)
    return chromosome, counts


def parallel_count_surrounding_SNPs(snp_df, flank, threads):
    processing_df = snp_df.copy()
    processing_df.sort_values(by=['chr', 'pos'], inplace=True)
    data_to_process = [(chromosome,
                        processing_df.loc[processing_df['chr'] == chromosome, 'pos'].values,
                        flank)
                       for chromosome in processing_df['chr'].unique()]
    pool = Pool(processes=threads)
    results = pool.map(count_surrounding_SNPs_chromosome, data_to_process)
    pool.close()
    pool.join()

    for chromosome, counts in results:
        processing_df.loc[processing_df['chr'] == chromosome, 'count'] = counts

    return processing_df

def filter_df_by_threshold(df, col_name, min_count, max_count):
    filtered_df = df[(df[col_name] >= min_count) & (df[col_name] <= max_count)]

    return filtered_df

if __name__ == '__main__':
    import argparse

    usage = """python Clustered_SNPs_filter.py -s SNPs_df.txt -f 50 -t 10 -o output"""

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument("-s", "--snp_df", dest="snp_df", action="store", nargs='?',
                        help="Input accessible SNPs dataframe (txt).", metavar="FILE", required=True)
    parser.add_argument("-f", "--flank", dest="flank", action="store", nargs='?', default="100",
                        help="Flanking length to filter out clustered SNPs (default: 100bp).",
                        metavar="STRING")
    parser.add_argument("--max", dest="max_snps", action="store", nargs='?', default="1",
                        help="Maximum number in flanking regions (default: 1).",
                        metavar="STRING")
    parser.add_argument("--min", dest="min_snps", action="store", nargs='?', default="1",
                        help="Minimum number in flanking regions (default: 1).",
                        metavar="STRING")
    parser.add_argument("-t", "--threads", dest="threads", action="store", nargs='?', default="1",
                        help="Threads (default: 1). The number of chromosomes is recommended.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()
    # 1 read SNPs dataframe
    print("Reading SNPs dataframe ...", datetime.datetime.now())
    headers = ["chr", "pos", "type", "ref", "alt"]
    SNPs_df = pd.read_csv(args.snp_df, delimiter='\t', names=headers, header=None)

    # 2 filter dataframe based on surrounding SNPs
    print("Filtering SNPs dataframe ...", datetime.datetime.now())
    count_df = parallel_count_surrounding_SNPs(SNPs_df, int(args.flank), int(args.threads))
    filtered_SNPs_df = filter_df_by_threshold(count_df, 'count', int(args.min_snps), int(args.max_snps))

    # 3 output SNPs dataframe
    output_SNPs_df = filtered_SNPs_df.drop('count', axis=1)
    output_SNPs_df.to_csv(f"{args.output}.cluster_filtered_SNPs.tsv", sep='\t', header=True, index=False)
    before_filtering = str(int(subprocess.check_output(["wc", "-l", args.snp_df]).split()[0]))
    after_filtering = str(diff_file_line_count(args.snp_df, f"{args.output}.cluster_filtered_SNPs.tsv"))
    print("Filtered out %s SNPs from %s SNPs in cluster filtering." % (after_filtering, before_filtering), datetime.datetime.now())
