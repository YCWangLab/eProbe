import datetime
from fasta_operator import read_fasta, write_fasta
from Get_probe_from_seq import get_probes_with_seq

def read_BED(BED_file):
    BED_dict = {}
    with open(BED_file) as intervals:
        for interval in intervals.readlines():
            interval = interval.strip()
            if interval:
                chrom, start, end = interval.split()[:3]
                if chrom not in BED_dict.keys():
                    BED_dict[chrom] = []
                    BED_dict[chrom].append([int(start), int(end)])
                else:
                    BED_dict[chrom].append([int(start), int(end)])
    for chrom, BED_list in BED_dict.items():
        BED_dict[chrom] = sorted(BED_list, key=lambda x: x[0])
    return BED_dict


def subtract_genome(fasta, BED):
    seq_dict = {}
    for chrom, blocks in BED.items():
        for block in blocks:
            seq_dict['_'.join([chrom, block[0], block[1]])] = fasta[chrom][int(block[0]) - 1:int(block[1])]
    return seq_dict

if __name__ == '__main__':
    import argparse

    usage = '''python Get_probe_from_bed.py -b *.bed -r ref.fasta  -l 100 -s 30 -o output'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-b', '--bed', dest='bed', action='store', nargs='?',
                        help='Input BED file for subtract sequences.',
                        metavar='FILE', required=True)
    parser.add_argument('-r', '--ref', dest='ref', action='store', nargs='?',
                        help='Input reference genome for subtracting sequences with BED.',
                        metavar='FILE', required=True)
    parser.add_argument('-l', '--len', dest='len', action='store', nargs='?',
                        help='Length of probes.',
                        metavar='STRING', required=True)
    parser.add_argument('-s', '--step', dest='step', action='store', nargs='?',
                        help='Step for generating probes.',
                        metavar='STRING', required=True)
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    genome = read_fasta(args.ref)
    BED = read_fasta(args.bed)
    templates = subtract_genome(genome, BED)

    output_probe_seq = {}
    for template_id, template in templates.items():
        probes = get_probes_with_seq(template, args.len, args.step)
        for i in range(1, len(probes) + 1, 1):
            output_probe_seq[f'{template_id}.{str(i)}'] = probes[i - 1]

    print("Generating probe sequences completed.", datetime.datetime.now())
    write_fasta(output_probe_seq, f'{args.output}.probe.fasta')