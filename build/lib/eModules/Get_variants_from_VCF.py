# -*- coding:utf-8 -*-
import vcf
import pandas as pd
import multiprocessing
import datetime


def read_genome_size(genome_size_file):
    genome_size_dict = {}
    with open(genome_size_file) as lines:
        for line in lines:
            chrom_id, chrom_len = line.strip().split("\t")
            genome_size_dict[chrom_id] = int(chrom_len)
    return genome_size_dict


def process_chrom_variants(chrom_vcf_reader, indel_mode, chrom):
    snp_dict = {"chr": [], "pos": [], "type": [], "ref": [], "alt": []}
    indel_sv_dict = {"chr": [], "start": [], "end": [], "type": []}
    for record in chrom_vcf_reader:
        if record.is_snp:
            pos = record.POS
            type_ = record.var_subtype
            ref = record.REF
            alt = record.ALT[0]
            snp_dict["chr"].append(chrom)
            snp_dict["pos"].append(pos)
            snp_dict["type"].append(type_)
            snp_dict["ref"].append(ref)
            snp_dict["alt"].append(alt)

        if indel_mode == "on" and (record.is_indel or record.is_sv):
            start = record.affected_start
            end = record.affected_end
            type_ = record.var_subtype
            indel_sv_dict["chr"].append(chrom)
            indel_sv_dict["start"].append(start)
            indel_sv_dict["end"].append(end)
            indel_sv_dict["type"].append(type_)

    snp_df = pd.DataFrame(snp_dict)
    indel_sv_df = pd.DataFrame(indel_sv_dict)
    return snp_df, indel_sv_df


def launch_process_chrom_variants(chrom):
    global vcf_file
    global indel_opt
    vcf_reader = vcf.Reader(filename=vcf_file)
    indel_opt = args.indel

    snp_df, indel_sv_df = process_chrom_variants(vcf_reader.fetch(chrom), indel_opt, chrom)
    return snp_df, indel_sv_df


if __name__ == "__main__":
    import argparse

    usage = """python Get_variants_from_VCF.py -g genome_size.txt -v *.vcf.gz --indel on -t threads -o output"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-g", "--g_size", dest="g_size", action="store", nargs='?',
                        help="Genome size file (tsv).",
                        metavar="FILE", required=True)
    parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                        help="Gzipped and indexed VCF file (vcf.gz).",
                        metavar="FILE", required=True)
    parser.add_argument("--indel_filter", dest="indel_filter", action="store", nargs='?',
                        help="Filter out those SNPs around structural variants in the VCF (on/off, default: on).",
                        choices=["on", "off"], default="on", metavar="STRING")
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default="1",
                        help="Number of thread (default: 1).", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()

    vcf_file = args.vcf_file

    print("Reading genome size file ...", datetime.datetime.now())
    genome_size = read_genome_size(args.g_size)
    chr_list = list(genome_size.keys())
    print("Processing VCF file in parallel...", datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = pool.map(launch_process_chrom_variants, [chrom for chrom in chr_list])
    pool.close()
    pool.join()

    print("Combining and outputting SNP and InDel tables...", datetime.datetime.now())
    Snp_df_processes, InDel_SV_df_processes = zip(*results)
    Snp_df_combined = pd.concat(Snp_df_processes, ignore_index=True)
    Snp_df_combined.sort_values(by=["chr", "pos"]).to_csv(f"{args.output}.processing_SNPs.tsv", "\t", index=False,
                                                          header=False)
    if args.indel == "on":
        InDel_SV_df_combined = pd.concat(InDel_SV_df_processes, ignore_index=True)
        InDel_SV_df_combined.sort_values(by=["chr", "start"]).to_csv(f"{args.output}.processing_InDels.tsv", "\t", index=False,
                                                                     header=False)
