from Bio.SeqUtils import MeltingTemp as Tm
from Bio.Seq import Seq
import pandas as pd
import multiprocessing
import datetime

from fasta_operator import split_fasta_into_dicts

def Tm_cal(seq, method, table):
    if method == 'tm_NN':
        if table == 'R_DNA_NN1':
            rna_seq = Seq(seq).transcribe()
            tm = Tm.Tm_NN(rna_seq, nn_table=Tm.R_DNA_NN1)
        elif table == 'DNA_NN1':
            tm = Tm.Tm_NN(seq, nn_table=Tm.DNA_NN1)
        elif table == 'DNA_NN2':
            tm = Tm.Tm_NN(seq, nn_table=Tm.DNA_NN2)
        elif table == 'DNA_NN3':
            tm = Tm.Tm_NN(seq, nn_table=Tm.DNA_NN3)
        elif table == 'DNA_NN4':
            tm = Tm.Tm_NN(seq, nn_table=Tm.DNA_NN4)
        else:
            raise ValueError('Chose the table for tm_NN from: DNA_NN1/DNA_NN2/DNA_NN3/DNA_NN4/R_DNA_NN1.')
    elif method == 'tm_GC':
        tm = Tm.Tm_GC(seq)
    else:
        raise ValueError('Chose the method from tm_NN or tm_GC for estimating melting temperature.')
    return tm


def process_probe_dict(probe_dict, method, table):
    tag_list = []
    for header, seq in probe_dict.items():
        chr = header.split(':')[0]
        start, end = header.split(':')[1].split('_')[0].split('-')
        pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
        tag = str(round(float(Tm_cal(seq, method, table)), 4))
        tag_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return tag_list


if __name__ == '__main__':
    import argparse

    usage = '''python Get_tag_tm.py -f probes.fasta -t threads -o output'''

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input fasta.', metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                        help='Number of threads (default: 1).',
                        default=1, metavar='STRING')
    parser.add_argument('-m', '--method', dest='method', action='store', nargs='?',
                        help='Method of calculating melting temperature (tm_NN/tm_GC; default: tm_NN).',
                        choices=['tm_NN', 'tm_GC'], default='tm_NN', metavar='STRING')
    parser.add_argument('--table', dest='table', action='store', nargs='?',
                        help='Table of thermodynamic NN values.'
                             'DNA/DNA hybridizations:'
                             'DNA_NN1: values from Breslauer et al. (1986)'
                             'DNA_NN2: values from Sugimoto et al. (1996)'
                             'DNA_NN3: values from Allawi & SantaLucia (1997)'
                             'DNA_NN4: values from SantaLucia & Hicks (2004)'
                             'For RNA/DNA hybridizations:'
                             'R_DNA_NN1: values from Sugimoto et al. (1995) (default).',
                        choices=['DNA_NN1', 'DNA_NN2', 'DNA_NN3', 'DNA_NN4', 'R_DNA_NN1'],
                        default='R_DNA_NN1', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)

    args = parser.parse_args()

    split_dicts = split_fasta_into_dicts(args.fasta, int(args.thread))
    print('Calculating melting temperature ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = [pool.apply_async(process_probe_dict, args=(fa_dict, args.method, args.table,)) for fa_dict in
                   split_dicts]
        all_SNPs_tag_list = []
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list, columns=['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt', 'tm'])
    SNPs_df.to_csv(args.output, sep='\t', index=False)
    print('Calculating melting temperature completed.', datetime.datetime.now())
