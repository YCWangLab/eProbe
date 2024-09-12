from Bio import SeqIO
import datetime
from fasta_operator import write_fasta

def get_probes_with_seq(seq, length, step):
    probes_list = []

    if len(seq) > length:
        for i in range(0, len(seq) - length + 1, step):
            probes_list.append(seq[i:i + length])
        if len(seq) % step != 0:
            probes_list.append(seq[-length:])
    return probes_list


if __name__ == '__main__':
    import argparse

    usage = '''python Get_probes_from_seq.py -f *.fasta -l length -s step -o output.fasta'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input sequences (*.fasta) for generating probes.',
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

    genes = list(SeqIO.parse(args.fasta, 'fasta'))

    probe_dict = {}

    for gene in genes:
        description = '\t'.join(gene.description.split()[1:]) if gene.description else ''
        probe_dict[gene.id] = {'description': description, 'seq': []}
        probe_dict[gene.id]['seq'] = get_probes_with_seq(str(gene.seq), int(args.len), int(args.step))

    output_probe_dict = {}
    for gene_id, probes in probe_dict.items():
        description = probes['description']
        for i in range(1, len(probes['seq']) + 1):
            header = f'{gene_id}.{i}' + ('\t' + description if description else '')
            output_probe_dict[header] = probes['seq'][i - 1]

    print("Generating probe sequences completed.", datetime.datetime.now())
    write_fasta(output_probe_dict, f'{args.output}.probe.fasta')


