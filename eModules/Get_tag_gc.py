from Bio.SeqUtils import gc_fraction
import pandas as pd
import multiprocessing
import datetime

from fasta_operator import split_fasta_into_dicts

def process_probe_dict(probe_dict):
    tag_list = []
    for header, seq in probe_dict.items():
        chr = header.split(':')[0]
        start, end = header.split(':')[1].split('_')[0].split('-')
        pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
        tag = str(round(float(gc_fraction(seq))*100, 4))
        tag_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return tag_list

if __name__ == '__main__':
    import argparse

    usage = ''' python Get_tag_gc.py -f probes.fasta -t threads -o output '''

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.', metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?', default=1,
                        help='Number of threads (default: 1).', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)

    args = parser.parse_args()

    split_dicts = split_fasta_into_dicts(args.fasta, int(args.thread))
    print('Calculating GC content ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = [pool.apply_async(process_probe_dict, args=(fa_dict,)) for fa_dict in split_dicts]
        all_SNPs_tag_list = []
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list, columns=['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt', 'gc'])
    SNPs_df.to_csv(args.output, sep='\t', index=False)
    print('Calculating GC content completed.', datetime.datetime.now())

