import subprocess
import pysam
import datetime
import os
from fasta_operator import read_fasta


def run_bowtie2(reference, k_, thread, mode, fasta, trim5, trim3, outpre):
    bowtie2_cmd = f'bowtie2 -f -x {reference} -k {k_} -p {thread} --{mode}-local --trim5 {trim5} --trim3 {trim3} -U {fasta} -S {outpre}.eProbe_temp.sam'
    try:
        subprocess.run(bowtie2_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print('Error running Bowtie2:', e)


def run_samtools(outpre):
    samtools_cmd = f'samtools view -F 4 -F 256 -Sb {outpre}.eProbe_temp.sam | samtools sort -o {outpre}.eProbe_temp.bam && samtools index {outpre}.eProbe_temp.bam'
    try:
        subprocess.run(samtools_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print('Error running samtools:', e)


def bam_filter_mode(bam_file, mode):
    mode = mode.lower()
    pass_ids = []

    with pysam.AlignmentFile(bam_file, 'rb') as bamfile:
        for alignment in bamfile:
            if alignment.is_supplementary:
                continue
            read_id = alignment.query_name
            mapq = alignment.mapping_quality
            nm_value = alignment.get_tag('NM')
            xo_value = alignment.get_tag('XO')

            if mode == 'strict':
                if (not alignment.has_tag('XS') and nm_value <= 3 and
                        mapq >= 30 and xo_value == 0):
                    pass_ids.append(read_id)
            elif mode == 'moderate':
                if nm_value <= 3 and mapq >= 30:
                    if not alignment.has_tag('XS'):
                        pass_ids.append(read_id)
                    else:
                        as_value = alignment.get_tag('AS')
                        xs_value = alignment.get_tag('XS')
                        if xs_value < as_value * 0.8:
                            pass_ids.append(read_id)

    return pass_ids


def bam_filter_custom(bam_file, para_dict):
    pass_ids = []

    min_mapq = int(para_dict.get('MAPQ', 30))
    max_nm = int(para_dict.get('NM', float(3)))
    max_xo = int(para_dict.get('XO', float(0)))
    max_xg = int(para_dict.get('XG', float(0)))
    max_xs_ratio = float(para_dict.get('XS_ratio', 0.8))

    with pysam.AlignmentFile(bam_file, 'rb') as bamfile:
        for alignment in bamfile:
            if alignment.is_supplementary:
                continue

            read_id = alignment.query_name
            mapq = alignment.mapping_quality
            nm_value = alignment.get_tag('NM')
            xo_value = alignment.get_tag('XO')
            xg_value = alignment.get_tag('XG')
            if not alignment.has_tag('XS'):
                if (mapq >= min_mapq and nm_value <= max_nm and xo_value <= max_xo and
                        xg_value <= max_xg):
                    pass_ids.append(read_id)
            else:
                as_value = alignment.get_tag('AS')
                xs_value = alignment.get_tag('XS')
                if (mapq >= min_mapq and nm_value <= max_nm and xo_value <= max_xo and
                        xg_value <= max_xg and xs_value < as_value * max_xs_ratio):
                    pass_ids.append(read_id)
    return pass_ids


def fasta_filter(fasta_dict, id_list, keep=True):
    filtered_dict = {}

    for seq_id, seq in fasta_dict.items():
        if (seq_id in id_list) == keep:
            filtered_dict[seq_id] = seq

    return filtered_dict


if __name__ == '__main__':
    import argparse

    usage = '''python Snp_filter_ac.py -f probe.fasta -i bowtie2_index -o output'''

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.',
                        metavar='FILE', required=True)
    parser.add_argument('-i', '--index', dest='index', action='store', nargs='?',
                        help='Bowtie2 database.',
                        metavar='FILE', required=True)
    parser.add_argument('-k', dest='keep_hit', action='store', nargs='?',
                        help='The number of hits to be kept in bowtie2 mapping (default: 100).',
                        default='100', metavar='STRING')
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                        help='Thread for running bowtie2 mapping.',
                        default='1', metavar='STRING')
    parser.add_argument('--mm', dest='m_mode', action='store', nargs='?',
                        help='Bowtie2 local mapping mode (very-fast, fast, sensitive, very-sensitive). '
                             '(default: sensitive)',
                        default='sensitive', metavar='STRING')
    parser.add_argument('--trim5', dest='trim5', action='store', nargs='?',
                        help='Trim bases from 5\'/left end of reads (default: 0).',
                        default='0', metavar='STRING')
    parser.add_argument('--trim3', dest='trim3', action='store', nargs='?',
                        help='Trim bases from 3\'/right end of reads (default: 0).',
                        default='0', metavar='STRING')
    parser.add_argument('--fm', dest='f_mode', action='store', nargs='?',
                        help='Filtering mode (strict, moderate, custom) (default: moderate).'
                             'If custom is enabled, custom filtering parameters will be used.',
                        default='moderate', metavar='STRING')
    parser.add_argument('--mapq', dest='mapq', action='store', nargs='?',
                        help='[custom mode] Minimum mapping quality (default: 30).',
                        default='30', metavar='STRING')
    parser.add_argument('--nm', dest='nm', action='store', nargs='?',
                        help='[custom mode] Maximum number of edit distances (default: 3).',
                        default='3', metavar='STRING')
    parser.add_argument('--xo', dest='xo', action='store', nargs='?',
                        help='[custom mode] Maximum gap number (default: 0).',
                        default='0', metavar='STRING')
    parser.add_argument('--xg', dest='xg', action='store', nargs='?',
                        help='[custom mode] Maximum total gap length (default: 0).',
                        default='0', metavar='STRING')
    parser.add_argument('--xsr', dest='xsr', action='store', nargs='?',
                        help='[custom mode] Maximum secondary alignment score ratio (default: 0.8).',
                        default='0.8', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='FILE', required=True)

    args = parser.parse_args()

    print('Launching bowtie2 ...', datetime.datetime.now())
    run_bowtie2(args.index, args.keep_hit, args.thread, args.m_mode.lower(), args.fasta, args.trim5, args.trim3,
                args.output)
    print('Mapping completed ...', datetime.datetime.now())

    print('Filtering SAM ...', datetime.datetime.now())
    run_samtools(args.output)
    print('Filtering completed ...', datetime.datetime.now())

    print('Start filtering SNPs in inaccessible regions...', datetime.datetime.now())
    if args.f_mode.lower() in ['strict', 'moderate']:
        pass_ids = bam_filter_mode(args.output + '.eProbe_temp.bam', args.f_mode.lower())
    elif args.f_mode.lower() == 'custom':
        parameters = {
            'MAPQ': args.mapq,
            'NM': args.nm,
            'XO': args.xo,
            'XG': args.xg,
            'XS_ratio': args.xsr,
        }
        pass_ids = bam_filter_custom(args.output + '.eProbe_temp.bam', parameters)
    else:
        raise ValueError('Invalid filtering mode specified.')

    probe_dict = read_fasta(args.fasta)
    filtered_probe_dict = fasta_filter(probe_dict, pass_ids, keep=True)
    with open(f'{args.output}.processing_probe.fasta', 'w') as out:
        for id, seq in filtered_probe_dict.items():
            out.write('>' + id + '\n')
            out.write(seq + '\n')
    out.close()

    rm_list = [f'{args.output}.eProbe_temp.sam', f'{args.output}.eProbe_temp.bam', f'{args.output}.eProbe_temp.bam.bai']
    for file in rm_list:
        os.remove(file)
