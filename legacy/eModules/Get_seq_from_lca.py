import argparse
from Bio import SeqIO
import re


def extract_modified_lines(file_path, patterns):
    modified_lines = []
    with open(file_path, 'r') as file:
        for line in file:
            if any(re.compile(r'\s+%s:' % pattern).search(line) for pattern in patterns.split(',')):
                parts = line.split('\t')[0].split(':')
                modified_line = ':'.join(parts[:-3])
                modified_lines.append(modified_line.strip())
    return set(modified_lines)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Extract sequences from a FastQ or Fasta file based on ngsLCA results(*.lca).
    Example: python Get_seq_from_lca.py -l *.lca -t 4530,4531,4532 -i *.fastq -f fasta  -o output.fasta ''')
    parser.add_argument('-l', '--lca', dest='lca', action='store', nargs='?',
                        help='ngsLCA result (*.lca)', metavar='FILE', required=True)
    parser.add_argument('-t', '--taxa', dest='taxa', action='store', nargs='?',
                        help='taxid to search for in the input lca (e.g., 4530 or 4530,4531,4532).',
                        metavar='STRING', required=True)
    parser.add_argument('-i', '--input_seq', dest='seq', action='store', nargs='?',
                        help='fastq or fastq format sequence file.', metavar='FILE', required=True)
    parser.add_argument('-f', '--format', dest='format', action='store', nargs='?', choices=['fasta', 'fastq'],
                        help='specify the seq file format.', required=True)
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?', metavar='STRING', help='Path to the output file.',
                        required=True)
    args = parser.parse_args()

    target_ids = extract_modified_lines(args.lca, args.taxa)

    with open(args.seq, 'r') as input_file:
        with open(args.output, 'w') as output_file:
            file_format = args.format.lower()
            for record in SeqIO.parse(input_file, file_format):
                if record.id in target_ids:
                    output_file.write(f'>{record.id}\n{record.seq.strip()}\n')
