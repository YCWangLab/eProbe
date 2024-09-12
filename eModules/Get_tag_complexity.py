import multiprocessing
import pandas as pd
import datetime

from fasta_operator import split_fasta_into_dicts

def cal_dust_score(seq):
    k = 3
    kmers_list = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmers_list.append(kmer)
    kmers_set = list(set(kmers_list))
    Ct_dict = {}
    for i in kmers_set:
        Ct_dict[i] = kmers_list.count(i)
    all_ct = 0
    for ct in Ct_dict.values():
        all_ct = all_ct + int(ct) * (int(ct) - 1) / 2
    return round(all_ct / (len(seq) - 3), 4)


def process_probe_dict(probe_dict):
    tag_list = []
    for header, seq in probe_dict.items():
        chr = header.split(':')[0]
        start, end = header.split(':')[1].split('_')[0].split('-')
        pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
        tag = str(cal_dust_score(seq))
        tag_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return tag_list


if __name__ == '__main__':

    import argparse

    usage = '''python Get_tag_complexity.py -f probes.fasta -t threads -o output'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.', metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?', default=1,
                        help='Number of threads (default: 1).', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)
    args = parser.parse_args()

    split_dicts = split_fasta_into_dicts(args.fasta, int(args.thread))
    print('Calculating complexity score ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = [pool.apply_async(process_probe_dict, args=(fa_dict,)) for fa_dict in split_dicts]
        all_SNPs_tag_list = []
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list,
                           columns=['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt', 'complexity'])
    SNPs_df.to_csv(args.output, sep='\t', index=False)
    print('Calculating complexity score completed.', datetime.datetime.now())
