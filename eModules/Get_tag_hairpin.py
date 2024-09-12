from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import pandas as pd
import collections
import multiprocessing
import datetime

from fasta_operator import split_fasta_into_dicts

def hairpin_PairwiseAligner_score(probe_dict, match=1, mismatch=-1, gap=-5, extend_gap=-5):
    hairpin_score = collections.OrderedDict()
    aligner = PairwiseAligner(mode='local', score='blastn')
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.gap_score = gap
    aligner.extend_gap_score = extend_gap
    for header, seq in probe_dict.items():
        alignments = aligner.align(seq, Seq(seq).reverse_complement())
        if len(alignments):
            best_alignment = max(alignments, key=lambda x: (x.score, len(x.target)))
            alignment_sim = str(best_alignment).split('\n')[1].strip()
            gaps = alignment_sim.count('-')
            matches = alignment_sim.count('|')
            mismatches = alignment_sim.count('.')
            # calculate score based on the number of matches, mismatches, and gaps
            score = 3 * int(matches) - 2 * int(mismatches) - 5 * int(gaps)
            hairpin_score[header] = score
    return hairpin_score


def process_probe_dict(probe_dict, match_s, mismatch_s, gap_s, extend_gap_s):
    tag_list = []
    hairpin = hairpin_PairwiseAligner_score(probe_dict, match_s, mismatch_s, gap_s, extend_gap_s)
    for header, score in hairpin.items():
        chr = header.split(':')[0]
        start, end = header.split(':')[1].split('_')[0].split('-')
        pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
        tag = str(score)
        tag_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return tag_list


if __name__ == '__main__':
    import argparse

    usage = '''python Get_tag_hairpin.py -f probes.fasta -t 16 -o output'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.', metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?', default=1,
                        help='Number of threads (default: 1)', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)
    parser.add_argument('--mp', dest='match_p', action='store', nargs='?',
                        help='Match score for PairwiseAlign (default: 1).',
                        default=1, metavar='STRING')
    parser.add_argument('--mmp', dest='mismatch_p', action='store', nargs='?',
                        help='Mismatch penalty for PairwiseAlign (default: -1).',
                        default=-1, metavar='STRING')
    parser.add_argument('--gp', dest='gap_p', action='store', nargs='?',
                        help='Gap open penalty for PairwiseAlign (default: -5).',
                        default=-5, metavar='STRING')
    parser.add_argument('--gpe', dest='extend_gap_p', action='store', nargs='?',
                        help='Gap extension  penalty for PairwiseAlign (default: -5).',
                        default=-5, metavar='STRING')
    args = parser.parse_args()

    split_dicts = split_fasta_into_dicts(args.fasta, int(args.thread))
    print('Calculating hairpin score ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = [pool.apply_async(process_probe_dict,
                                    args=(fa_dict, int(args.match_p), int(args.mismatch_p), int(args.gap_p),
                                          int(args.extend_gap_p),)) for fa_dict in split_dicts]
        all_SNPs_tag_list = []
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list,
                           columns=['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt', 'hairpin'])
    SNPs_df.to_csv(args.output, sep='\t', index=False)
    print('Calculating hairpin score completed.', datetime.datetime.now())
