# -*- coding:utf-8 -*-
import vcf
import pandas as pd
import argparse
from Bio import SeqIO
import collections

usage = """python Get_Accessible_SNPs.py -f genome.fasta -v *.vcf.gz -d 100 -b genome_blocks.bed -o output.txt"""


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-f", "--genome_fasta", dest="genome_fasta", action="store", nargs='?',
                    help="Reference genome fasta", metavar="FILE")
parser.add_argument("-v", "--raw_vcf", dest="vcf_file", action="store", nargs='?',
                    help="Input filtered vcf.gz file, with tbi index", metavar="FILE")
parser.add_argument("-d", "--distance", dest="distance", action="store", nargs='?',
                    help="Distance for filtering SNPs around InDels/SVs", metavar="STRING")
parser.add_argument("-b", "--bed", dest="bed", action="store", nargs='?', default=None,
                    help="Genome regions that need to retain their SNPs (bed format)", metavar="FILE")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs table", metavar="FILE")
args = parser.parse_args()

## parameters
# genome fasta
fasta_file = args.genome_fasta
# filtered vcf file
vcf_file = args.vcf_file
# distance around InDels/SVs
dist = int(args.distance)
# bed format
access_file = args.bed
# output SNPs table
outsnp = args.output


# function for reading genome
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict


# function for creating a snp_df
def make_snp_df(vcf_reader_, chr_list_):
    snp_dict = {"chr": [], "pos": [], "type": []}
    for chrom in chr_list_:
        for record in vcf_reader_.fetch(chrom):
            if record.is_snp:
                # one-based coordinate
                pos = record.POS
                # ts, tv
                type = record.var_subtype
                snp_dict["chr"].append(chrom)
                snp_dict["pos"].append(pos)
                snp_dict["type"].append(type)
    snp_df = pd.DataFrame(snp_dict)
    return snp_df


# function for creating an indel_df (bed format)
def make_InDels_SVs_df(vcf_reader_, chr_list_):
    InDels_SVs_dict = {"chr": [], "start": [], "end": [], "type": []}
    for chrom in chr_list_:
        for record in vcf_reader_.fetch(chrom):
            if record.is_indel or record.is_sv:
                # one-based coordinate
                start = record.affected_start
                end = record.affected_end
                # ins, del
                type = record.var_subtype
                InDels_SVs_dict["chr"].append(chrom)
                InDels_SVs_dict["start"].append(start)
                InDels_SVs_dict["end"].append(end)
                InDels_SVs_dict["type"].append(type)
    InDels_SVs_df = pd.DataFrame(InDels_SVs_dict)
    return InDels_SVs_df


# function for removing snp around blocks
def remove_snp_by_bed(snp_df_, bed_df_, distance):
    snp_drop_list = []
    for row in bed_df_.itertuples():
        chrom = row.chr
        start = int(row.start) - int(distance)
        end = int(row.end) + int(distance)
        # index list of snp around indel and sv
        drop_list = snp_df_[(snp_df_.chr == chrom) & (snp_df_.pos > start) & (snp_df_.pos < end)].index.tolist()
        if len(drop_list):
            snp_drop_list = snp_drop_list + drop_list
    snp_drop_list = sorted(list(set(snp_drop_list)))
    drop_snp_df = snp_df_.drop(snp_drop_list)
    sort_snp_df = drop_snp_df.sort_index()
    return sort_snp_df


# function for keeping snp around blocks
def keep_snp_by_bed(snp_df_, bed_df_, distance):
    snp_keep_list = []
    for row in bed_df_.itertuples():
        chrom = row.chr
        start = int(row.start) - int(distance)
        end = int(row.end) + int(distance)
        # index list of snp within blocks
        keep_list = snp_df_[(snp_df_.chr == chrom) & (snp_df_.pos > start) & (snp_df_.pos < end)].index.tolist()
        if len(keep_list):
            snp_keep_list = snp_keep_list + keep_list
    snp_keep_list = sorted(list(set(snp_keep_list)))
    return snp_df_.loc[snp_keep_list]


# function for making normal bed df
def make_normal_bed_df(bed_file_):
    # one-based coordinate
    names_col = ["chr", "start", "end"]
    bed_df = pd.read_table(bed_file_, sep="\t", names=names_col,
                           dtype={"start": int, "end": int})
    return bed_df


# read genome
genome = read_fasta(fasta_file)
# make chr list
chr_list = list(genome.keys())
# read vcf
vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
# make snp df
snp_df_raw = make_snp_df(vcf_reader, chr_list)
# make indel df
indel_sv_df = make_InDels_SVs_df(vcf_reader, chr_list)
# remove snp around indel
filtered_snp_df = remove_snp_by_bed(snp_df_raw, indel_sv_df, dist)
if access_file:
    # read accessible region file
    accessible_df = make_normal_bed_df(access_file)
    # keep snp within accessible regions
    snp_df_acc = keep_snp_by_bed(filtered_snp_df, accessible_df, 0)
else:
    snp_df_acc = filtered_snp_df

# snp_df_acc contains all accessible SNPs in both coding regions and non-coding regions
snp_df_acc.to_csv(outsnp, "\t", index=False)
# output genome length information
with open("Chr_length.txt", "w") as out:
    for k, v in genome.items():
        out.write(k + "\t" + str(len(v)) + "\n")
out.close()
