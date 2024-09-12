import pandas as pd
import collections
import datetime
import random
from fasta_operator import read_fasta


def get_probes_with_df(snp_df, method, fa_dict, left_flank=0, right_flank=0, shift=0):
    probe_dict = collections.OrderedDict()
    for row in snp_df.itertuples():
        chromosome = row.chr
        if method == 'pos':
            # 1-based to 0-based
            start = int(row.pos) - left_flank - 1 + int(shift)
            end = int(row.pos) + right_flank - 1 + int(shift)
            ref_start = max(0, start)
            ref_end = min(end, len(fa_dict[chromosome]) - 1)
            header = f'{chromosome}:{int(ref_start + 1)}-{int(ref_end + 1)}_{row.pos}_{row.type}_{row.ref}_{row.alt}'
        elif method == 'edge':
            # 1-based to 0-based
            ref_start = int(row.start) - 1
            ref_end = int(row.end) - 1
            header = f'{chromosome}:{int(row.start)}-{int(row.end)}_{row.pos}_{row.type}_{row.ref}_{row.alt}'
        else:
            raise ValueError('Unsupported method specified. Use "pos" or "edge".')
        sequence = fa_dict[chromosome][int(ref_start):int(ref_end) + 1]
        probe_dict.setdefault(chromosome, {})[header] = sequence.upper()

    return probe_dict


def replace_targetSNP(probe_dict):
    processed_probe_dict = {}
    base_set = {'A', 'T', 'C', 'G'}
    for chr, chr_seq in probe_dict.items():
        processed_probe_dict[chr] = {}
        for header, seq in chr_seq.items():
            start, end = header.split(':')[1].split('_')[0].split('-')
            pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
            replace_pos = int(float(pos)) - int(float(start))
            new_base = random.choice(list(base_set - {ref, alt}))
            replaced_probe = seq[:replace_pos] + new_base + seq[replace_pos + 1:]
            processed_probe_dict[chr][header] = replaced_probe.upper()

    return processed_probe_dict


def filter_type(probe_dict, select_type):
    processed_probe_dict = {}
    for chr, chr_seq in probe_dict.items():
        processed_probe_dict[chr] = {}
        for header, seq in chr_seq.items():
            pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
            if type_ == select_type:
                processed_probe_dict[chr][header] = seq
    return processed_probe_dict


if __name__ == '__main__':
    import argparse

    usage = '''python Get_probe_from_snp.py -d SNPs_df.txt -g genome.fasta -l 81 -s shift -o output'''

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument('-d', '--snp_df', dest='snp_df', action='store', nargs='?',
                        help='Input SNPs dataframe (tsv).', metavar='FILE', required=True)
    parser.add_argument('-g', '--genome', dest='genome', action='store', nargs='?',
                        help='Reference genome file (fasta).', metavar='FILE', required=True)
    parser.add_argument('-l', '--length', dest='length', action='store', nargs='?', default=81,
                        help='Length of probes (odd number recommended).', metavar='STRING')
    parser.add_argument('-s', '--shift', dest='shift', action='store', nargs='?', default='0',
                        help='Position shift (default: 0) fo method pos. - for shifting left and + for shifting right.',
                        metavar='STRING')
    parser.add_argument('-m', '--method', dest='method', action='store', nargs='?',
                        help='Method [pos/edge] for extracting seq from genome (default: pos).',
                        choices=['pos', 'edge'], default='pos', metavar='STRING')
    parser.add_argument('--replace', dest='replace', action='store', nargs='?',
                        help='Replace target SNP with a base neither ref not alt (on/off, default: off).',
                        choices=['on', 'off'], default='off', metavar='STRING')
    parser.add_argument('--type', dest='type', action='store', nargs='?',
                        help='Select probes based on mutation type (ts/tv, default: both).',
                        choices=['ts', 'tv'], default=None, metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='FILE', required=True)

    args = parser.parse_args()

    print('Reading SNPs dataframe ...', datetime.datetime.now())
    SNPs_df = pd.read_csv(args.snp_df, delimiter='\t', header=0, index_col=False)

    print('Reading reference genome: %s ...' % args.genome, datetime.datetime.now())
    genome = read_fasta(args.genome)

    print('Extracting probe sequences ...', datetime.datetime.now())
    if int(args.length) % 2 == 1:
        left = right = int(args.length) // 2
    else:
        left = (int(args.length) / 2) - 1
        right = left + 1
    if abs(int(args.shift)) < min(right, left):
        probes_dict = get_probes_with_df(SNPs_df, args.method, genome, left, right, int(args.shift))
    else:
        raise ValueError('The target SNPs were shifted out of probes. Please reduce the length of your shift.')

    if args.replace == 'on':
        print('Replacing target SNPs ...', datetime.datetime.now())
        probes_dict = replace_targetSNP(probes_dict)

    if args.type:
        print(f'Select {args.type} SNPs ...', datetime.datetime.now())
        probes_dict = filter_type(probes_dict, args.type)

    probe_num = sum(len(chr_probe) for chr_probe in probes_dict.values())
    print(f'Generated {probe_num} probes.', datetime.datetime.now())

    if args.method == 'pos' and args.replace == 'off':
        output_path = f'{args.output}.processing_probe.fasta'
    else:
        output_path = f'{args.output}.filtered_probe.fasta'
    with open(output_path, 'w') as out_file:
        for chr, chr_seq in probes_dict.items():
            for header, seq in chr_seq.items():
                out_file.write(f'>{header}\n{seq}\n')
