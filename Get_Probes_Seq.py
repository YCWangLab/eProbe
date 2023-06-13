# -*- coding:utf-8 -*-
import pandas as pd
import argparse
import collections
from Bio import SeqIO

usage = """python Get_Probes_Seq.py -s SNPs_df.txt -g genome.fasta -l 81 -o output.fasta"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-d", "--snp_df", dest="snp_df", action="store", nargs='?',
                    help="Input accessible SNPs dataframe (.txt)", metavar="FILE")
parser.add_argument("-g", "--genome", dest="genome_fasta", action="store", nargs='?',
                    help="Reference genome file in fasta format", metavar="FILE")
parser.add_argument("-l", "--length", dest="probes_length", action="store", nargs='?',
                    help="Length of probes (odd number recommanded)", metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output fasta file storing probe sequences", metavar="FILE")

args = parser.parse_args()

SNPs_df_file = args.snp_df
genome_file = args.genome_fasta
Probe_len = int(args.probes_length)
output_file = args.output


# function: read reference fasta
# note: vcf_pos is 1-based, while ref_pos (in python) is 0-based, so vcf_pos = ref_pos + 1
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict


# function: extract flanks based on a given SNPs df
def get_flanks_with_df(snp_df_, left, right, fa_dict):
    probes_dict = collections.OrderedDict()
    for row in snp_df_.itertuples():
        # 1-based
        chr = row.chr
        if chr not in probes_dict.keys():
            probes_dict[chr] = {}
        pos = row.pos
        # to 0-based
        if int(pos) - left - 1 > 0:
            ref_start = int(pos) - left - 1
        else:
            ref_start = 0
        if int(pos) + right - 1 < len(fa_dict[chr]):
            ref_end = int(pos) + right - 1
        else:
            ref_end = len(fa_dict[chr])
        # make header
        header = chr + ":" + str(ref_start) + "-" + str(ref_end) + "_" + str(row.pos) + "_" + row.type
        # extract sequence
        seq = fa_dict[chr][ref_start:(ref_end + 1)]
        probes_dict[chr][header] = seq
    return probes_dict


# 1 read SNPs DF
SNPs_df = pd.read_csv(SNPs_df_file, delimiter='\t', header=0)

# 2 read genome
genome = read_fasta(genome_file)

# 3 get probes
# for odd number (recommanded)
if Probe_len % 2 == 1:
    left = int((Probe_len - 1) / 2)
    right = left
# for even number
elif Probe_len % 2 != 1:
    left = int(Probe_len / 2) - 1
    right = left + 1
probes_dict = get_flanks_with_df(SNPs_df, left, right, genome)
# 4 output probes
with open(output_file, "w") as out:
    for chr, chr_seq in probes_dict.items():
        for id, seq in chr_seq.items():
            out.write(">" + id + "\n")
            out.write(seq + "\n")
out.close()
