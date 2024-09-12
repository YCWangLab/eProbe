import pysam
import pandas as pd
import multiprocessing
import datetime
import warnings
from fasta_operator import read_chr_size


def is_snp(record):
    if len(record.ref) != 1:
        return False
    for alt in record.alts:
        if alt not in ['A', 'C', 'G', 'T']:
            return False
    return True


def is_transition(ref, alt):
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    return (ref, alt) in transitions or (alt, ref) in transitions


def determine_tv_ts(record):
    ref = record.ref
    alt = record.alts[0]

    if is_transition(ref, alt):
        return 'ts'
    else:
        return 'tv'


def process_chrom_snps(chrom_vcf_path, chrom):
    snp_dict = {'chr': [], 'pos': [], 'type': [], 'ref': [], 'alt': []}
    non_biallelic = 0
    with pysam.VariantFile(chrom_vcf_path) as vcf:
        for record in vcf.fetch(chrom):
            if is_snp(record):
                if len(record.alts) == 1:
                    snp_dict['chr'].append(chrom)
                    snp_dict['pos'].append(record.pos)
                    snp_dict['type'].append(determine_tv_ts(record))
                    snp_dict['ref'].append(record.ref)
                    snp_dict['alt'].append(record.alts[0])
                else:
                    non_biallelic = non_biallelic + 1
    if non_biallelic:
        warnings.warn(f'{non_biallelic} non-biallelic SNPs detected in chromosome {chrom}.', Warning)
    else:
        print(f'No non-biallelic SNPs detected in chromosome {chrom}.')

    return pd.DataFrame(snp_dict)


def launch_process_chrom_snps(args_tuple):
    chrom, vcf_path = args_tuple
    return process_chrom_snps(vcf_path, chrom)


if __name__ == '__main__':
    import argparse

    usage = '''python Get_snp_from_vcf.py -g chromosome_size.tsv -v *.vcf.gz -t threads -o output'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                        help='Gzipped and indexed VCF file (vcf.gz).',
                        metavar='FILE', required=True)
    parser.add_argument('-g', '--genome_size', dest='genome_size', action='store', nargs='?',
                        help='Chromosome size file (tsv).',
                        metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                        help='Number of thread (default: 1).',
                        default='1', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    print('Reading chromosome size file ...', datetime.datetime.now())
    chr_size = read_chr_size(args.genome_size)
    chr_list = list(chr_size.keys())
    print('Retrieving SNPs from VCF. ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = pool.map(launch_process_chrom_snps, [(chrom, args.vcf) for chrom in chr_list])

    pool.close()
    pool.join()

    snp_df_combined = pd.concat(results, ignore_index=True)
    print(f'Retrieved {len(snp_df_combined)} SNPs from VCF.', datetime.datetime.now())
    snp_df_combined.sort_values(by=['chr', 'pos']).to_csv(f'{args.output}.processing_SNPs.tsv', '\t', index=False,
                                                          header=False)
