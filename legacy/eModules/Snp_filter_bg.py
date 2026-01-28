import subprocess
import os
import datetime
from pipe_controler import diff_file_line_count


def run_kraken2(fasta, databse, thread, min_hit, confidence, unclassified_fasta):
    kraken2_cmd = f'kraken2 --db {databse} --threads {thread} --minimum-hit-groups {min_hit} --unclassified-out {unclassified_fasta} --output - --confidence {confidence} {fasta}'
    subprocess.run(kraken2_cmd, shell=True, check=True)


if __name__ == '__main__':
    import argparse

    usage = '''python Snp_filter_bg.py -f probe.fasta -d kraken2_databse(s) 
    -t 10 --min_hits 2 -o output'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.',
                        metavar='FILE', required=True)
    parser.add_argument('-d', '--db', dest='db', action='store', nargs='?',
                        help='Kraken2 databases (separated with comma) for background filtering.',
                        metavar='STRING', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                        help='Number of thread (default: 1).',
                        default=1, metavar='STRING')
    parser.add_argument('--min_hits', dest='min_hits', action='store', nargs='?',
                        help='Minimum number of non-overlapping kmer hits (separate kmer stacks) '
                             'to filter out (default: 2).',
                        default=2, metavar='STRING')
    parser.add_argument('--confidence', dest='confidence', action='store', nargs='?',
                        help='Minimum percentage of kmers out of a sequence needed to be assigned '
                             'to a node (default: 0).',
                        default=0, metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    print('Start background filtering ...', datetime.datetime.now())
    output_file = f'{args.output}.processing_probe.fasta'
    if len(args.db.split(',')) == 1:
        run_kraken2(args.fasta, args.db, int(args.thread), args.min_hits, args.confidence,
                    output_file)
    else:
        temp_input = args.fasta
        temp_output = f'{args.output}.background_filter_tempout'
        for db in args.db.split(','):
            run_kraken2(temp_input, db, int(args.thread), args.min_hits, args.confidence,
                        temp_output)
            temp_input = f'{args.output}.background_filter_tempin'
            os.rename(temp_output, temp_input)
        os.rename(temp_input, output_file)
    before_filtering = str(int(int(subprocess.check_output(['wc', '-l', args.fasta]).split()[0]) / 2))
    after_filtering = str(int(diff_file_line_count(args.fasta, output_file) / 2))
    print('Filtered out %s from %s probes in background noise filtering.' % (after_filtering, before_filtering),
          datetime.datetime.now())
