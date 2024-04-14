from Bio.Align import PairwiseAligner
from Bio import SeqIO
import pandas as pd
import collections
import multiprocessing
import os
import shutil
from fasta_operator import split_fasta,read_fasta

def hairpin_PairwiseAligner_score(probes, match=1, mismatch=-1, gap=-5, extend_gap=-5):
    hairpin_score = collections.OrderedDict()
    aligner = PairwiseAligner(mode='local', score="blastn")
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.gap_score = gap
    aligner.extend_gap_score = extend_gap
    for record in probes:
        alignments = aligner.align(record.seq, record.seq, strand="-")
        if len(alignments):
            # get best alignment
            best_alignment = max(alignments, key=lambda x: (x.score, len(x.target)))
            # parse the alignment
            alignment_sim = str(best_alignment).split("\n")[1].strip()
            gaps = alignment_sim.count("-")
            matches = alignment_sim.count("|")
            mismatches = alignment_sim.count(".")
            # calculate score based on the number of matches, mismatches and gaps
            score = 3 * int(matches) - 2 * int(mismatches) - 5 * int(gaps)
            hairpin_score[record.id] = score
    return hairpin_score


def process_probe_file(probe_file, match_s, mismatch_s, gap_s, extend_gap_s):
    probes_record = list(SeqIO.parse(probe_file, "fasta"))
    hairpin = hairpin_PairwiseAligner_score(probes_record, match_s, mismatch_s, gap_s, extend_gap_s)
    score_list = []
    for header, score in hairpin.items():
        chr = header.split(":")[0]
        start, end = header.split(":")[1].split("_")[0].split("-")
        pos, type_, ref, alt = header.split(":")[1].split("_")[1:5]
        tag = str(score)
        score_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return score_list


if __name__ == '__main__':
    import argparse

    usage = """python Get_tag_hairpin.py -f probes.fasta -t 16 -o output"""

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Input fasta.", metavar="FILE", required=True)
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of threads (default: 1)", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)
    parser.add_argument("--mp", dest="match_p", action="store", nargs='?', default=1,
                        help="Match score for PairwiseAlign (default: 1).",
                        metavar="STRING")
    parser.add_argument("--mmp", dest="mismatch_p", action="store", nargs='?', default=-1,
                        help="Mismatch penalty for PairwiseAlign (default: -1).", metavar="STRING")
    parser.add_argument("--gp", dest="gap_p", action="store", nargs='?', default=-5,
                        help="Gap open penalty for PairwiseAlign (default: -5).", metavar="STRING")
    parser.add_argument("--gpe", dest="extend_gap_p", action="store", nargs='?', default=-5,
                        help="Gap extension  penalty for PairwiseAlign (default: -5).", metavar="STRING")
    args = parser.parse_args()

    split_files = split_fasta(args.fasta, int(args.thread), args.output)

    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        all_SNPs_tag_list = []
        results = [pool.apply_async(process_probe_file,
                                    args=(probe_file, args.match_p, args.mismatch_p, args.gap_p, args.extend_gap_p,))
                   for probe_file in split_files]
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list,
                           columns=["chr", "start", "end", "pos", "type", "ref", "alt", "hairpin"])
    SNPs_df.to_csv(args.output + ".tsv", sep='\t', index=False)

    shutil.rmtree(args.output + "_eProbe_temp.hairpin")
