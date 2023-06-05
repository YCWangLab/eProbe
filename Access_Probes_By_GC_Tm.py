# -*- coding:utf-8 -*-
import collections
import pandas as pd
import argparse
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as Tm

usage = """python Assess_Probes_By_GC_Tm.py -s SNPs_df.txt -g genome.fasta -k 71,121,10 -m Tm_NN"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-d", "--snp_df", dest="snp_df", action="store", nargs='?',
                    help="Input accessible SNPs dataframe (.txt)", metavar="FILE")
parser.add_argument("-g", "--genome", dest="genome_fasta", action="store", nargs='?',
                    help="Reference genome file in fasta format", metavar="FILE")
parser.add_argument("-k", "--kmers", dest="kmers", action="store", nargs='?',
                    help="Kmers for generating probes. e.g., 71,121,10 (start,end,step) or 81 (odd number)", metavar="STRING")
parser.add_argument("-m", "--method", dest="Tm_method", action="store", nargs='?',
                    help="Method for melting temperature calculation (Tm_GC/Tm_NN)", metavar="STRING")

args = parser.parse_args()

SNPs_df_file = args.snp_df
genome_file = args.genome_fasta
input_kmer = args.kmers
if len(input_kmer.split(",")) == 3:
    k_start = int(input_kmer.split(",")[0])
    k_end = int(input_kmer.split(",")[1])
    k_step = int(input_kmer.split(",")[2])
    kmer_list = [i for i in range(k_start, k_end + k_step, k_step)]
elif len(input_kmer.split(",")) == 1:
    kmer_list = []
    kmer_list.append(int(input_kmer))
flank_lens = [math.ceil((int(i) - 1) / 2) for i in kmer_list]
Tm_method = args.Tm_method

# function: read reference fasta
# note: vcf_pos is 1-based, while ref_pos (in python) is 0-based, so vcf_pos = ref_pos + 1
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict


# function: extract flanks based on a given SNPs df
def get_flanks_with_df(snp_df_, flank_len, fa_dict):
    probes_dict = collections.OrderedDict()
    for row in snp_df_.itertuples():
        # 1-based
        chr = row.chr
        if chr not in probes_dict.keys():
            probes_dict[chr] = {}
        pos = row.pos
        # to 0-based
        if int(pos) - int(flank_len) - 1 > 0:
            ref_start = int(pos) - int(flank_len) - 1
        else:
            ref_start = 0
        if int(pos) + int(flank_len) - 1 < len(fa_dict[chr]):
            ref_end = int(pos) + int(flank_len) - 1
        else:
            ref_end = len(fa_dict[chr])
        # make header
        header = chr + ":" + str(ref_start) + "-" + str(ref_end)
        # extract sequence
        seq = fa_dict[chr][ref_start:(ref_end + 1)]
        probes_dict[chr][header] = seq
    return probes_dict


# reading SNPs dataframe
SNPs_df = pd.read_csv(SNPs_df_file, delimiter='\t', header=0)
# reading genome fasta
genome = read_fasta(genome_file)
# making GC and Tm distribution
for flank in flank_lens:
    probes = get_flanks_with_df(SNPs_df, flank, genome)
    GC_list = []
    Tm_list = []
    for chr_probes in probes.values():
        for probe in chr_probes.values():
            gc = round(GC(probe), 4)
            if Tm_method == "Tm_NN":
                #  Allawi & SantaLucia (1997)
                tm = Tm.Tm_NN(probe, nn_table=Tm.DNA_NN3)
            elif Tm_method == "Tm_GC":
                tm = Tm.Tm_GC(probe)
            GC_list.append(gc)
            Tm_list.append(tm)
    GC_mean = np.mean(GC_list)
    GC_std = np.std(GC_list)
    Tm_mean = np.mean(Tm_list)
    Tm_std = np.std(Tm_list)
    sns.distplot(GC_list, color="green")
    sns.distplot(Tm_list)
    plt.title("Distribution of GC (green) and Tm (Blue)")
    plt.savefig('Distribution_of_GC_and_Tm_for_K%s.pdf' % str(flank * 2 + 1), dpi=800)
    plt.close()
    print("K = %s" % str(flank * 2 + 1))
    print("GC mean: %s" % str(GC_mean))
    print("GC std: %s" % str(GC_std))
    print("Tm mean: %s" % str(Tm_mean))
    print("Tm std: %s" % str(Tm_std))
