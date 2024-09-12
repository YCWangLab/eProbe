from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import os
import datetime
from fasta_operator import read_fasta, write_fasta


def multiple_sequence_alignment(sequences, temp_base):
    input_file = f'{temp_base}.temp.fasta'
    output_file = f'{temp_base}.temp.aln'
    temp_file = f'{temp_base}.temp.dnd'
    write_fasta(sequences, input_file)

    # run clustalw
    clustalw_cline = ClustalwCommandline('clustalw2', infile=input_file)
    stdout, stderr = clustalw_cline()

    if os.path.exists(output_file):
        alignment = AlignIO.read(output_file, 'clustal')
    else:
        raise FileNotFoundError(f'ClustalW failed to produce output file: {output_file}')

    # remove temp
    os.remove(input_file)
    os.remove(output_file)
    os.remove(temp_file)

    return alignment


def sliding_window_variation(alignment, window_size=100, step_size=30):
    alignment_length = alignment.get_alignment_length()
    variations = []

    for start in range(0, alignment_length - window_size + 1, step_size):
        end = start + window_size
        window = alignment[:, start:end]
        variation_count = 0
        for i in range(end - start):
            column = window[:, i]
            if len(set(column)) > 1:  # check difference (mutation)
                variation_count += 1

        variations.append((start, end, variation_count, window))

    return variations


def find_original_position(aligned_seq, original_seq):
    # remove gap '-' and upper the seq
    trimmed_seq = aligned_seq.replace('-', '').upper()

    # upper the seq before searching
    original_position = original_seq.upper().find(trimmed_seq)

    return original_position


def get_probes_with_seq(seq, length, step):
    probes_list = []
    i = 0
    if len(seq) > length:
        while True:
            if i + int(length) < len(seq):
                probe = seq[i:i + int(length)]
                if len(probe) == length:
                    probes_list.append(probe)
                i = i + int(step)
            else:
                probe = seq[-int(length):]
                if len(probe) == length:
                    probes_list.append(probe)
                break
    return probes_list


def get_probes_with_alleles(alleles_dict, temp_base, window_size, step_size):
    probe_seq = {}

    for gene, alleles in alleles_dict.items():
        if gene not in probe_seq.keys():
            probe_seq[gene] = []

        allele_seqs = {}
        for num, allele_seq in enumerate(alleles.values()):
            allele_seqs[str(num)] = allele_seq
        # have alleles
        if len(allele_seqs) > 1:
            # align different alleles
            alignment = multiple_sequence_alignment(allele_seqs, temp_base)
            # detect variation from alignment
            window_alignments = sliding_window_variation(alignment, window_size, step_size)
            for start, end, variation_count, window_alignment in window_alignments:
                # no variation
                if int(variation_count) == 0:
                    probe_seq[gene].append(str(window_alignment[0].seq))
                # have variation
                else:
                    alleles_in_window = []
                    for record in window_alignment:
                        start = int(find_original_position(str(record.seq), str(allele_seqs[record.id])))
                        end = start + window_size
                        alleles_in_window.append(str(allele_seqs[record.id])[start:end])
                    probe_seq[gene] = probe_seq[gene] + list(set(alleles_in_window))

        else:
            probe_seq[gene] = get_probes_with_seq(list(allele_seqs.values())[0], window_size, step_size)

    return probe_seq


def process_fasta(fasta_dict, method, separator, min_freq):
    alleles_dict = {}

    def add_to_dict(gene_id, allele_num, seq):
        if gene_id not in alleles_dict:
            alleles_dict[gene_id] = {}
        alleles_dict[gene_id][allele_num] = seq

    for header, seq in fasta_dict.items():
        if method.lower() == 'bed':
            gene_id = header.split('_Allele')[0]
            allele_num = header.split('_Allele')[1].split('_Freq')[0]
            freq = header.split('_Freq')[1]
            if float(freq) >= float(min_freq):
                add_to_dict(gene_id, allele_num, seq)
        elif method.lower() == 'seq':
            gene_id, allele_num = header.split(separator)[:2]
            add_to_dict(gene_id, allele_num, seq)

    return alleles_dict


if __name__ == '__main__':
    import argparse

    usage = '''python Get_probe_from_allele.py -f *.fasta -l 100 -s 30 -o output.fasta'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('--input', dest='input', action='store', nargs='?',
                        help='Source of input data (seq/bed, default: bed).'
                             'seq: allele sequences collected by user ("seq" Method in Seq generator).'
                             'bed: allele sequences inferred by eProbe ("bed" Method in Seq generator).',
                        choices=['seq', 'bed'], default='bed', metavar='STRING')
    parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                        help='Input template sequences (*.fasta) for generating probes.',
                        metavar='FILE', required=True)
    parser.add_argument('-l', '--len', dest='len', action='store', nargs='?',
                        help='Length of probes.',
                        metavar='STRING', required=True)
    parser.add_argument('-s', '--step', dest='step', action='store', nargs='?',
                        help='Step for generating probes.',
                        metavar='STRING', required=True)
    parser.add_argument('--min_freq', dest='min_freq', action='store', nargs='?',
                        help='Minimum frequency for alleles to be considered (default: 0.05).'
                             'Only use when --input is "bed".',
                        default=0.05, metavar='STRING')
    parser.add_argument('--sep', dest='sep', action='store', nargs='?',
                        help='Separator for indicating allele num (default: _). '
                             'Only use when --input is "seq".',
                        default='_', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    fasta_dict = read_fasta(args.fasta)

    alleles_dict = process_fasta(fasta_dict, args.input, args.sep, float(args.min_freq))
    print('Generating probe sequences ...', datetime.datetime.now())
    probe_seq = get_probes_with_alleles(alleles_dict, 'eProbe', int(args.len), int(args.step))

    output_probe_seq = {}
    for gene, probes in probe_seq.items():
        for num, probe in enumerate(probes, 1):
            seq_id = f'{gene}.{num}'
            output_probe_seq[seq_id] = probe

    print('Generating probe sequences completed.', datetime.datetime.now())
    write_fasta(output_probe_seq, f'{args.output}.probe.fasta')
