from Bio.SeqUtils import GC
import pandas as pd
import collections
from Bio import SeqIO
import multiprocessing
from fasta_operator import read_fasta

def process_probe(header, seq, func):
    chr = header.split(":")[0]
    start, end = header.split(":")[1].split("_")[0].split("-")
    pos, type_, ref, alt = header.split(":")[1].split("_")[1:5]
    tag = str(round(func(seq), 4))
    return [chr, start, end, pos, type_, ref, alt, tag]


if __name__ == '__main__':
    import argparse

    usage = """ python Get_tag_GC.py -f probes.fasta -t threads -o output """

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Input fasta.", metavar="FILE", required=True)
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of threads (default: 1).", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()

    probes = read_fasta(args.fasta)
    SNPs_tag_list = []
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = []
        for header, seq in probes.items():
            results.append(pool.apply_async(process_probe, (header, seq, GC)))
        for result in results:
            SNPs_tag_list.append(result.get())
        sorted_SNPs_tag_list = sorted(SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list, columns=["chr", "start", "end", "pos", "type", "ref", "alt", "GC"])
    SNPs_df.to_csv(args.output + ".tsv", sep='\t', index=False)
