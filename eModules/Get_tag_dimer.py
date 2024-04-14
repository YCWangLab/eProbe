from Bio import SeqIO
import os
import collections
import multiprocessing
import datetime
import pandas as pd
import shutil
from fasta_operator import split_fasta

def reverse_complement(fasta):
    rcseq_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta, "fasta"))
    for record in records:
        rcseq_dict[record.id] = str(record.seq.replace("\n", "").reverse_complement())
    return rcseq_dict


def kmer_processor_seq(seq, k, kmer_freq_dict):
    n = len(seq)
    for i in range(n - k + 1):
        kmer = seq[i:i + k]
        kmer_freq_dict[kmer] = kmer_freq_dict.get(kmer, 0) + 1


def kmer_processor_dict(fasta_dict, k):
    kmer_freq_dict = {}
    for seq in fasta_dict.values():
        kmer_processor_seq(seq, k, kmer_freq_dict)
    return kmer_freq_dict


def merge_kmer_dicts(*dicts):
    merged_dict = {}
    for dictionary in dicts:
        for key, value in dictionary.items():
            if key in merged_dict:
                merged_dict[key] += value
            else:
                merged_dict[key] = value
    return merged_dict


def filter_dict_by_value(dictionary, min_value):
    df = pd.DataFrame.from_dict(dictionary, orient='index').reset_index()
    df.columns = ['key', 'value']

    df = df[df['value'] >= int(min_value)]

    filtered_dict = df.set_index('key')['value'].to_dict()

    return filtered_dict


def cal_dimer_score_seq(seq, k, kmer_freq_dict):
    seq_kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    dimer_score = 0

    for kmer in seq_kmers:
        if kmer in kmer_freq_dict:
            dimer_score += kmer_freq_dict[kmer]
    return dimer_score


def cal_dimer_score_fasta(fasta_file, k, kmer_freq_dict):
    dimer_score = {}
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    for record in sequences:
        id = record.id
        score = cal_dimer_score_seq(str(record.seq), k, kmer_freq_dict)
        dimer_score[id] = score
    return dimer_score


def process_probe_file(probe_file, k, kmer_freq_dict):
    dimer_score = cal_dimer_score_fasta(probe_file, k, kmer_freq_dict)
    result_list = []
    for header, score in dimer_score.items():
        chr = header.split(":")[0]
        start, end = header.split(":")[1].split("_")[0].split("-")
        pos, type_, ref, alt = header.split(":")[1].split("_")[1:5]
        tag = str(score)
        result_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return result_list


if __name__ == '__main__':
    import argparse

    usage = """python Get_tag_dimer.py -f fasta -k kmer_length --min_freq 2 -t thread -o output_prefix"""

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Fasta file for count kmer frequency.", metavar="FILE", required=True)
    parser.add_argument("-k", "--kmer", dest="kmer", action="store", nargs='?', default=11,
                        help="Length of kmer (default: 11).", metavar="STRING")
    parser.add_argument("--min_freq", dest="min_freq", action="store", nargs='?', default=2,
                        help="Minimum value of frequency.", metavar="STRING")
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of thread (default: 1).", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)
    args = parser.parse_args()

    print("Splitting fasta ...", datetime.datetime.now())
    split_files = split_fasta(args.fasta, int(args.thread) * 2, args.output)

    print("Counting reverse complement kmer frequency ...", datetime.datetime.now())
    rc_fasta_list = []
    for split_file in split_files:
        rc_fasta_list.append(reverse_complement(split_file))
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = []
        for split_rc_fasta in rc_fasta_list:
            results.append(pool.apply_async(kmer_processor_dict, args=(split_rc_fasta, int(args.kmer))))
        rc_kmer_freq_dict = {}
        for result in results:
            rc_kmer_freq_dict = merge_kmer_dicts(rc_kmer_freq_dict, result.get())

    print("Filtering reverse complement kmer with frequency lower than %s ..." % str(args.min_freq),
          datetime.datetime.now())
    filtered_kmer_dict = filter_dict_by_value(rc_kmer_freq_dict, args.min_freq)

    print("Calculating dimer score ...", datetime.datetime.now())
    all_SNPs_tag_list = []
    results = []
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        for split_file in split_files:
            results.append(pool.apply_async(process_probe_file, args=(split_file, int(args.kmer), filtered_kmer_dict)))
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list,
                           columns=["chr", "start", "end", "pos", "type", "ref", "alt", "dimer"])
    SNPs_df.to_csv(args.output + ".tsv", sep='\t', index=False)
    print("Calculating dimer score finished.", datetime.datetime.now())

    if os.path.exists(args.output + "_eProbe_temp.dimer"):
        shutil.rmtree(args.output + "_eProbe_temp.dimer")
